
library(ipw)
library(ltmle)
library(dplyr)
library(tidyr) 
library(ggplot2)
library(survival)
library(DataCombine)

#######################
# Data Preparation
#######################

data(sampleDataForLtmleMSM)
dt <- sampleDataForLtmleMSM$data

data.long <- reshape(dt, direction = "long", 
                     varying = list(c("A0", "A1", "A2"), c("CD4_0", "CD4_1", "CD4_2"), c("Y1", "Y2", "Y3")), 
                     timevar = "time", times = 0:2, 
                     v.names = c("A", "CD4", "Y"), 
                     idvar = c("patientID"))
data.long <- arrange(data.long, data.long$patientID, data.long$time)

n <- max(data.long$patientID)
T <- max(data.long$time)

for (i in 0:T){
  if(i==0){
    data.long$Z[data.long$time == i] <- rbinom(n, size = 1, prob = 0.05)
    data.long$Y[data.long$time == i] <- (1-data.long$Z[data.long$time == i]) * data.long$Y[data.long$time == i]
    data.long$atrisk[data.long$time == i] <- 1
  }else{
    data.long$Z[data.long$time == i] <- (1-data.long$Y[data.long$time == i-1]) * pmax(data.long$Z[data.long$time == i-1], rbinom(n, size = 1, prob = 0.05))
    data.long$Y[data.long$time == i] <- (1-data.long$Z[data.long$time == i]) * data.long$Y[data.long$time == i]
    data.long$atrisk[data.long$time == i] <- 1 - pmax(data.long$Y[data.long$time == i-1], data.long$Z[data.long$time == i-1])
  }
}

# Apply LOCF for missing treatment indicators
data.long <- fill(group_by(data.long, patientID), A)


# Set CD4 of inactive patients to NA
data.long$CD4[data.long$atrisk == 0] <- NA

# Create lagged 1 variable for the treatment start indicator
data.long <- slide(data.long, Var = "A", GroupVar = "patientID", NewVar = "A_LAG1", slideBy = -1)
data.long$A_LAG1 <- ifelse(is.na(data.long$A_LAG1), 0, data.long$A_LAG1)

df <- data.long[,c(7, 3, 1, 2, 5, 4, 10, 6, 8, 9)]

### Data description
# patientID - patient indicator
# time - time
# male - patient gender
# age - patient age
# CD4 - patient CD4 count
# A - treatment indicator
# A_LAG1 - LAG1 treatment indicator
# Y - event of interest
# Z - competing event
# atrisk - indicator of being "at risk"

# Note: for time = 0, the corresponding Y and Z in that row, should be interpreted as Y_1 and Z_1, 
# and more generally for the row encoding time k, it should be interpreted as Y_{k+1} and Z_{k+1} 
# as is usually done when time k would correspond with tstart (of the time interval) and time k+1 
# with the end of the time interval (or the start of the next time interval)

################################################################################
# Define considered regime *thresholds*
################################################################################

CD4_regime <- c(330, 350, 370)

################################################################################
# Propensity score model
################################################################################

ps_model <- glm(formula = A ~ time + male + age + CD4,
                data = df[df$A_LAG1 == 0 & df$atrisk == 1,],
                family = "binomial")
summary(ps_model)

df$ps[df$atrisk == 1 & df$A_LAG1 == 0] <- predict(ps_model, df = df[df$atrisk == 1 & df$A_LAG1 == 0,], type = "response")
df$ps[df$atrisk == 1 & df$A_LAG1 == 1] <- 1
df$ps[df$atrisk == 0] <- NA 

################################################################################
# Create extended data set
################################################################################

# Create copy of a dataset for each regime
df_regimes <- data.frame(df[rep(seq_len(nrow(df)), times = length(CD4_regime)), ], CD4_regime = rep(CD4_regime, each = nrow(df)))

################################################################################
# Create treatment *counterfactual* initiation indicator based on different treatment regimes
################################################################################

# Treatment start indicator based on given treatment regime decision rule
df_regimes <- fill(mutate(group_by(df_regimes, patientID, CD4_regime), d = cummax(CD4 < CD4_regime)), d) #*** this makes sure the last observation for a patient is carried forward in case of missing values (e.g. if d = 1 before death or discharge, then d = 1 after death/discharge; same should hold for A)

# Factorize regime variable
df_regimes$CD4_regime <- as.factor(df_regimes$CD4_regime)

################################################################################
# Create compatibility indicator for different treatment regimes
################################################################################

# Compute compatibility indicator (comp) with given regime
df_regimes <- mutate(group_by(df_regimes, patientID, CD4_regime), comp = cumprod(A==d))

################################################################################
### Unstabilised IPCW
################################################################################

# Compute (unstabilized) Inverse Probability of Compatibility Weights using propensity scores and decision rules for different regimes. 
df_regimes <- mutate(group_by(df_regimes, patientID, CD4_regime), w_IPW = ifelse(comp == 0, 0, 1 / cumprod(ps^d * (1-ps)^(1-d)))) 

# Apply LOCF for missing IPC weights
df_regimes <- fill(group_by(df_regimes, patientID, CD4_regime), w_IPW)

# Make plots for obtained unstabilised IP weights (for compatible patients of a given regime)
ipwplot(df_regimes$w_IPW[df_regimes$comp == 1], timevar = df_regimes$time[df_regimes$comp == 1], logscale = FALSE, binwidth = 1, xlab = "Time", ylab = "IPCW")

################################################################################
### Aalen-Johanson estimator
################################################################################

df_regimes$T_start <- df_regimes$time
df_regimes$T_stop <- df_regimes$time + 1
df_regimes$event <- df_regimes$Y * 1 + df_regimes$Z * 2
df_regimes$event <- factor(df_regimes$event, 0:2, c("atrisk", "Y", "Z")) 

# Keep only compatible and at risk patients
df_regimes_comp_atrisk <- df_regimes[df_regimes$atrisk==1 & df_regimes$comp==1,] 

# Compute cumulative incidence curves based on AJ estimator using unstabilised weights
ci_1 <- survfit(Surv(T_start, T_stop, event) ~ 1, data = df_regimes_comp_atrisk[df_regimes_comp_atrisk$CD4_regime==CD4_regime[1],], weights = w_IPW, id = patientID)
ci_2 <- survfit(Surv(T_start, T_stop, event) ~ 1, data = df_regimes_comp_atrisk[df_regimes_comp_atrisk$CD4_regime==CD4_regime[2],], weights = w_IPW, id = patientID)
ci_3 <- survfit(Surv(T_start, T_stop, event) ~ 1, data = df_regimes_comp_atrisk[df_regimes_comp_atrisk$CD4_regime==CD4_regime[3],], weights = w_IPW, id = patientID)
ci_obs <- survfit(Surv(T_start, T_stop, event) ~ 1, data = df_regimes[df_regimes$CD4_regime==unique(df_regimes$CD4_regime)[1],], id = patientID)

################################################################################
### Marginal Structural Models (MSMs)
################################################################################

# MSM model with unstabilised IPCW for event of interest. 
# No conditioning on baseline covariates in the hazard model
MSM_Y <- glm(formula = Y ~ time * CD4_regime, 
             data = df_regimes_comp_atrisk[df_regimes_comp_atrisk$Z==0,],
             weights = df_regimes_comp_atrisk[df_regimes_comp_atrisk$Z==0,]$w_IPW,
             family = "binomial")

# Predicted ICU death hazards from the MSM model
df_regimes_comp_atrisk$Y_hazard_IPW <- predict(MSM_Y, df_regimes_comp_atrisk, type = "response")

################################################################################
### MSM models for competing event
################################################################################

# MSM model with unstabilised IPW for competing event. 
# No conditioning on baseline covariates in the hazard model
MSM_Z <- glm(formula = Z ~ time * CD4_regime, 
             data = df_regimes_comp_atrisk,
             weights = df_regimes_comp_atrisk$w_IPW,
             family = "binomial")

# Predicted ICU death hazards from the MSM model 
df_regimes_comp_atrisk$Z_hazard_IPW <- predict(MSM_Z, df_regimes_comp_atrisk, type = "response")

################################################################################
### Computation of Cumulative Incidence based on the formula (3) using MSM hazards 
################################################################################

# Computation based on the formula (3). Compute probabilities for each patient (based on his/her covariates)
df_regimes_comp_atrisk <- mutate(group_by(df_regimes_comp_atrisk, CD4_regime, patientID), prod = (1-Y_hazard_IPW)*(1-Z_hazard_IPW), S = cumprod((1-Y_hazard_IPW)*(1-Z_hazard_IPW))) 
df_regimes_comp_atrisk <- mutate(df_regimes_comp_atrisk, P = cumsum(Y_hazard_IPW*(1-Z_hazard_IPW)*lag(S, default = 1))) 

# Average across all patients per time wave and per regime 
df_MSM_res <- summarise(group_by(df_regimes_comp_atrisk, CD4_regime, time), P = mean(P))
df_MSM_res <- mutate(df_MSM_res, time = time + 1)
df_MSM_res.wide <- reshape(as.data.frame(df_MSM_res), direction = "wide", 
                           timevar = "CD4_regime", 
                           v.names = "P", 
                           idvar = "time")

# The estimator is 0 at time 0. (We estimate deaths for time 1 onwards).
df_MSM_res.wide <- rbind(rep(0, length(df_MSM_res.wide)), df_MSM_res.wide)

################################################################################
### Plot cumulative incidence curves for different treatment regimes
################################################################################

# Cumulative incidence curves based on the MSMs
plot(df_MSM_res.wide[,1], df_MSM_res.wide[,2], "s", xlab = "Time", ylab = "Cumulative incidence", col = "grey")
lines(df_MSM_res.wide[,1], df_MSM_res.wide[,3], "s", col = "blue")
lines(df_MSM_res.wide[,1], df_MSM_res.wide[,4], "s", col = "green")

# Cumulative incidence curves based on the Aalen-Johanson estimator
lines(ci_1[,2], fun = "event", conf.int = TRUE, col = "grey", lwd = 2)
lines(ci_2[,2], fun = "event", conf.int = TRUE, col = "blue", lwd = 2)
lines(ci_3[,2], fun = "event", conf.int = TRUE, col = "green", lwd = 2)

# Cumulative incidence curve for the observed cumulative incidence
lines(ci_obs[,2], fun = "event", conf.int = TRUE, col = "black", lwd = 2)

legend("bottomright", legend = c("CD4 < 330", "CD4 < 350", "CD4 < 370", "obs"), bty = "n", lty = 1, lwd = 2, col = c("grey", "blue", "green", "black"), cex = 0.9)

################################################################################
### Plot proportion of patients receiving treatment under particular DTR
################################################################################

A_number <- c()
d_number <- c()

# Select the data for one of the treatment regimes
data <- df_regimes[df_regimes$CD4_regime==CD4_regime[1],]
time <- unique(data$time)

# Calculate the proportion of patients receiving treatment at each time point under the current standard of care and under the considered DTR
for (i in time){
  A_number <- append(A_number, length(unique(data$patientID[data$A == 1 & data$time == i]))/length(unique(data$patientID[data$time == 0])))
  d_number <- append(d_number, sum(data$w_IPW[which(data$d == 1 & data$comp == 1 & data$time == i)])/sum(data$w_IPW[data$comp==1 & data$time == i]))
}

# Plot the proportion of patients receiving treatment at each time point under the current standard of care and under the considered DTR
data_plot <- cbind(0:2, A_number, d_number)
data_plot <- as.data.frame(data_plot)
names(data_plot) <- c("Time", "obs", "DTR")

data_plot_long <- reshape(data_plot, 
                          direction = "long",
                          v.names = "Value",
                          varying = c("obs", "DTR"),
                          times = c("obs", "DTR"))
names(data_plot_long) <- c("Time", "Regime", "Value", "id")
data_plot_long$Regime <- as.factor(data_plot_long$Regime)

ggplot(data_plot_long, aes(x=Time, y=Value, group=Regime, ymin = 0, ymax = 1)) +
  geom_point(aes(shape=Regime, color=Regime), size=3) +
  xlab("Time") + 
  ylab("Proportion of patients receiving treatment") 



