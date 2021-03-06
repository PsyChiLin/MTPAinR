######################### Permutation Test : 1D temporal clustering (pt_1dtc) #########################

### Load Supporting Tools
# Use pacman package to manage all the packages
library(pacman)
pacman::p_load(ERP, mnormt, fdrtool,
               tidyverse, gridExtra, crayon, 
               boot, reshape2, ggthemes, 
               devtools,randomForest,leaps, pROC)
# clear the working directory
rm(list = ls())

### Read LIFG Data
# Read the used Data, which is collected by the present NIRS experiment, in the present study
dtaAll <- readRDS("Data/NIRSdata_LTFGLMTG.Rdata")
# Define time point
tp <- seq(0,16,by=0.0959)
# Use the LIFG data only for demonstration purpose
dta <- filter(dtaAll, Area == "LIFG") 

### Explore Data
head(dta)
str(dta)
dim(dta)

### Start MUA_pt_1dtc

### Generate null distribution of suprathreshold cluster size (max STCS) 
permdta <- dta
# Resapmpling
permnum <- 10000
# Create a matrix to store Maximum Statistic: max STCS
permrst <- matrix(NA,permnum,1) 
# for loop for null distribution
for (j in 1:permnum){
  # Print to see the progress
  print(j)
  # Set seed for reproduciblility
  set.seed(j) 
  # Record the statistics at each sampled time point
  rdrst <- matrix(NA,2,170)
  # Randomly shuffle "Subject and Condition"
  permdta[,c(1:3)] <- permdta[sample(28),c(1,2,3)]
  # From 1 to 167 time points
  for (i in 4:170){
    # Each sampled time point
    testdta <- permdta[,c(1,2,3,i)]
    colnames(testdta)[4] <- "TimeN"
    # Fit a paired T model
    test <- t.test(TimeN~Condition, data = testdta, paired = T)
    # Store the t value and p value into the matrix
    rdrst[1,i] <- test$statistic
    rdrst[2,i] <- test$p.value
  }
  # Store the Maximum Statistic: max STCS into the matrix
  r <- rle(rdrst[2,4:170] < 0.05)
  if (TRUE %in% r$values){
    permrst[j,1] <- max(r$lengths[r$values == TRUE])
  } else {
    permrst[j,1] <-  0
  }
}

### Save the result to results folder as "MUA_pt_1dtc_Rst.Rdata"
# saveRDS(permrst,"Results/MUA_pt_1dtc_Rst.Rdata")
MUA_pt_1dtc_Rst <- readRDS("Results/MUA_pt_1dtc_Rst.Rdata")

### Show the max STCS and the significance threshold
#MUA_pt_1dtc_Rst[,1][is.na(MUA_pt_1dtc_Rst[,1])] <- 0
maxstat_STCS <- quantile(MUA_pt_1dtc_Rst[,1], probs = 0.95)
maxstat_STCS

### Results: Real data versus Significance threshold
# Read MUA results
MUA_Rst <- readRDS("Results/MUA_Rst.Rdata")

# Is there any significant time point ?
cluster <- rle(MUA_Rst$pvalue < 0.05) # or 0.025 in a setting of two tails
max(cluster$lengths[cluster$values == TRUE]) > maxstat_STCS

# The biggest cluster in the real dataset is 23(alpha = 0.025: better control in a two tails setting) or 37 (alpha = 0.05: less control in a one tail setting).
# Both not greater than maxstat_STCS, 37.
# max(cluster$lengths[cluster$values == TRUE])



