######################### Mass Univariate Analysis in GLM setting #########################

### Load Supporting Tools
# Use pacman package to manage all the packages
library(pacman)
pacman::p_load(ERP, mnormt, fdrtool,
               tidyverse, gridExtra, crayon, 
               boot, reshape2, ggthemes, 
               devtools,randomForest,leaps, pROC)
# clear the working directory
rm(list = ls())

### Read Data
# Read the used Data, which is collected by the present NIRS experiment, in the present study
dtaAll <- readRDS("Data/NIRSdata_LTFGLMTG.Rdata")
# Define time point
tp <- seq(0,16,by=0.0959)
# LIFG data
LIFG <- filter(dtaAll, Area == "LIFG")
# LMTG data
LMTG <- filter(dtaAll, Area == "LMTG")

### Start Mass Univariate Analysis in GLM setting (Similar to MTPA for comparisions)
### Set parameters of GLM
# Once a time point (MAU)
binwidth = 1
# Cross validation times
rcvnum <- 100
# Confidence interval
ci <- c(0.05,0.95)
# Set the uppr and lower bound 
upperbound <- 167-binwidth+1
lowerbound <- 1-binwidth+1
# Store the results
rst_LIFG <- matrix(NA,6,upperbound)
rst_LMTG <- matrix(NA,6,upperbound)

# LIFG : Start GLM model fitting
for (i in lowerbound:upperbound){
  # Print the progress
  if (i %in% c(seq(10,170,by =10))) {print(i)}
  # Record the AUC and CE
  ceauc <- matrix(NA,4,rcvnum)
  # Start cross validation 
  for (k in 1:rcvnum){
    # Set seed for reproduciable research
    set.seed(k+2018)
    # Training and Testing Data
    idc_test <- c(sample(1:14,5),sample(15:28,5))
    idc_train <- -idc_test
    # Fit an GLM model
    fit <- glm(Condition~.,data = LIFG[idc_train,c(3,(i+3):(i+3+binwidth-1))],
               family = binomial(link = "logit"))
    yhat_test <- predict(fit,newdata = LIFG[idc_test,],type = "response")
    yhat_test_class  <- ifelse(yhat_test < 0.5, 1, 0)
    ce_test <- mean(yhat_test_class!=LIFG[idc_test,]$Condition)
    auc_test <- pROC::auc(LIFG[idc_test,]$Condition,yhat_test)
    ceauc[2,k] <- ce_test
    ceauc[4,k] <- auc_test
  }
  # Store the results of CV
  rst_LIFG[1,i] <- mean(ceauc[2,])
  rst_LIFG[2,i] <- mean(ceauc[4,])
  rst_LIFG[3,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[1]
  rst_LIFG[4,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[2]
  rst_LIFG[5,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[1]
  rst_LIFG[6,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[2]
}
# Reorganize the results
LIFGm <- as.data.frame(rst_LIFG)
colnames(LIFGm) <- paste0("Time",1:length(1:upperbound))
row.names(LIFGm) <- c("CE","AUC","CE_l","CE_u","AUC_l","AUC_u")

# LMTG : Start GLM model fitting
for (i in lowerbound:upperbound){
  if (i %in% c(seq(10,170,by =10))) {print(i)}
  ceauc <- matrix(NA,4,rcvnum)
  for (k in 1:rcvnum){
    set.seed(k+2018)
    idc_test <- c(sample(1:14,5),sample(15:28,5))
    idc_train <- -idc_test
    fit <- glm(Condition~.,data = LMTG[idc_train,c(3,(i+3):(i+3+binwidth-1))],
               family = binomial(link = "logit"))
    yhat_test<- predict(fit,newdata = LMTG[idc_test,],type = "response")
    yhat_test_class  <- ifelse(yhat_test < 0.5, 1, 0)
    ce_test <- mean(yhat_test_class!=LMTG[idc_test,]$Condition)
    auc_test <- pROC::auc(LMTG[idc_test,]$Condition,yhat_test)
    ceauc[2,k] <- ce_test
    ceauc[4,k] <- auc_test
  }
  rst_LMTG[1,i] <- mean(ceauc[2,])
  rst_LMTG[2,i] <- mean(ceauc[4,])
  rst_LMTG[3,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[1]
  rst_LMTG[4,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[2]
  rst_LMTG[5,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[1]
  rst_LMTG[6,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[2]
}
LMTGm <- as.data.frame(rst_LMTG)
colnames(LMTGm) <- paste0("Time",1:length(1:upperbound))#
row.names(LMTGm) <- c("CE","AUC","CE_l","CE_u","AUC_l","AUC_u")

### Combine teh results of LIFG and LMTG
# LIFG
temp <- as.data.frame(t(LIFGm))
temp$Area <- "LIFG"
temp$Times <- tp
# LMTG
temp2 <- as.data.frame(t(LMTGm))
temp2$Area <- "LMTG"
temp2$Times <- tp
# Combine
MAU_GLM_Rst <- rbind(temp,temp2)
# Remove unused data
rm(temp,temp2,rst_LMTG,rst_LIFG,LIFG,LMTG,LIFGm,LMTGm,ceauc)

### Save the result to results folder as "MAU_GLM_Rst.Rdata"
# saveRDS(MAU_GLM_Rst,"Results/MAU_GLM_Rst.Rdata")
# MAU_GLM_Rst<- readRDS("Results/MAU_GLM_Rst.Rdata")

### Plot the results quickly
# ggplot(data = MAU_GLM_Rst, aes(x =Times, y = AUC, col = Area))+
#   geom_line()+
#   facet_grid(~Area)+
#   geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area, col = Area),
#               alpha = 0.3)+
#   theme_bw()





