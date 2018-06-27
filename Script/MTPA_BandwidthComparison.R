######################### Supplementary MTPA Binwidth Comparision #########################

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
# Read the bin2 and bin3 data
b2 <- readRDS("Results/MTPA_bin2_Rst.Rdata")
b3 <- readRDS("Results/MTPA_bin3_Rst.Rdata")

### Mean and SD
mean(b2$AUC[1:167]);sd(b2$AUC[1:167])
mean(b3$AUC[1:167]);sd(b3$AUC[1:167])

### T.test to compare bandwidth
t.test(b2$AUC[1:167],b3$AUC[1:167],var.equal = T)

# t = -0.65341, df = 332, p-value = 0.5139
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval: -0.019607022  0.009829377
# sample estimates:
#   mean of x mean of y 
# 0.7416198 0.7465086 