######################### Multi-Time Points Analysis (MTPA) #########################

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

### Start MTPA
### Set parameters for MTPA
# Consider 2 time points at each testing
binwidth = 2 
# Cross validation times
rcvnum <- 3
# Confidence interval
ci <- c(0.05,0.95)
# Set the upper and lower bound 
upperbound <- 167-binwidth+1
lowerbound <- 1
# Store the results
rst_LIFG <- matrix(NA,6,upperbound)
rst_LMTG <- matrix(NA,6,upperbound)

### LIFG : Start MTPA model fitting with RF
for (i in lowerbound:upperbound){
  # Print the progress
  if (i %in% c(seq(10,170,by =10))) {print(i)}
  # Record the AUC and CE
  ceauc <- matrix(NA,4,rcvnum)
  # Start cross validation
  for (k in 1:rcvnum){
    # Set seed for reproduciable research
    set.seed(k)
    # Training and Testing Data
    idc_test <- c(sample(1:14,5),sample(15:28,5))
    idc_train <- -idc_test
    # Fit an RF model
    fit <- randomForest(Condition~.,data = LIFG[idc_train,c(2,3,(i+3):(i+3+binwidth-1))],importance = F)
    yhat_test_prob <- predict(fit,newdata = LIFG[idc_test,],type = "prob")[,2]
    yhat_test_class <- predict(fit,newdata = LIFG[idc_test,],type = "class")
    # Record the results of RF fitting on Testing data
    ce_test <- mean(yhat_test_class!=LIFG[idc_test,]$Condition)
    auc_test <- pROC::auc(LIFG[idc_test,]$Condition,yhat_test_prob)
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
# Reorganize the results, avaerage all the time points that used to estimate the results
LIFGm <- matrix(NA,6,167)
LIFGm[,1] <- rst_LIFG[,1]
LIFGm[,167] <- rst_LIFG[,166]
for (i in 1:(upperbound-1)){
  tpi <- i + 1
  LIFGm[1,tpi] <- mean(rst_LIFG[1,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LIFGm[2,tpi] <- mean(rst_LIFG[2,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LIFGm[3,tpi] <- mean(rst_LIFG[3,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LIFGm[4,tpi] <- mean(rst_LIFG[4,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LIFGm[5,tpi] <- mean(rst_LIFG[5,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LIFGm[6,tpi] <- mean(rst_LIFG[6,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
}
LIFGm <- as.data.frame(LIFGm)
colnames(LIFGm) <- paste0("Time",1:167)
row.names(LIFGm) <- c("CE","AUC","CE_l","CE_u","AUC_l","AUC_u")

### LMTG : Start MTPA model fitting with RF
for (i in lowerbound:upperbound){
  # Print the progress
  if (i %in% c(seq(10,170,by =10))) {print(i)}
  # Record the AUC and CE
  ceauc <- matrix(NA,4,rcvnum)
  # Start cross validation
  for (k in 1:rcvnum){
    # Set seed for reproduciable research
    set.seed(k)
    # Training and Testing Data
    idc_test <- c(sample(1:14,5),sample(15:28,5))
    idc_train <- -idc_test
    # Fit an RF model
    fit <- randomForest(Condition~.,data = LMTG[idc_train,c(2,3,(i+3):(i+3+binwidth-1))],importance = F)
    yhat_test_prob <- predict(fit,newdata = LMTG[idc_test,],type = "prob")[,2]
    yhat_test_class <- predict(fit,newdata = LMTG[idc_test,],type = "class")
    # Record the results of RF fitting on Testing data
    ce_test <- mean(yhat_test_class!=LMTG[idc_test,]$Condition)
    auc_test <- pROC::auc(LMTG[idc_test,]$Condition,yhat_test_prob)
    ceauc[2,k] <- ce_test
    ceauc[4,k] <- auc_test
  }
  # Store the results of CV
  rst_LMTG[1,i] <- mean(ceauc[2,])
  rst_LMTG[2,i] <- mean(ceauc[4,])
  rst_LMTG[3,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[1]
  rst_LMTG[4,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[2]
  rst_LMTG[5,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[1]
  rst_LMTG[6,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[2]
}
# Reorganize the results, avaerage all the time points that used to estimate the results
LMTGm <- matrix(NA,6,167)
LMTGm[,1] <- rst_LMTG[,1]
LMTGm[,167] <- rst_LMTG[,166]
for (i in 1:(upperbound-1)){
  tpi <- i + 1
  LMTGm[1,tpi] <- mean(rst_LMTG[1,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LMTGm[2,tpi] <- mean(rst_LMTG[2,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LMTGm[3,tpi] <- mean(rst_LMTG[3,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LMTGm[4,tpi] <- mean(rst_LMTG[4,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LMTGm[5,tpi] <- mean(rst_LMTG[5,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
  LMTGm[6,tpi] <- mean(rst_LMTG[6,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
}
LMTGm <- as.data.frame(LMTGm)
colnames(LMTGm) <- paste0("Time",1:167)
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
MTPA_bin2_Rst <- rbind(temp,temp2)
# Remove unused data
rm(temp,temp2,rst_LMTG,rst_LIFG,LIFG,LMTG,LIFGm,LMTGm,ceauc)

### Save the result to results folder as "MTPA_bin2_Rst.Rdata"
# saveRDS(MTPA_bin2_Rst, "Results/MTPA_bin2_Rst.Rdata")
# MTPA_bin2_Rst <- readRDS("Results/MTPA_bin2_Rst.Rdata")

### Plot the results quickly
ggplot(data = MTPA_bin2_Rst,aes(x =Times, y = AUC, col = Area))+
  geom_line(size = 1.2)+
  geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area, col = Area),alpha = 0.3)+
  theme_bw()+
  facet_grid(~Area)
ggplot(data = MTPA_bin2_Rst,aes(x =Times, y = AUC, col = Area))+
  geom_line(size = 1.2)+
  geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area, col = Area),alpha = 0.3)+
  theme_bw()