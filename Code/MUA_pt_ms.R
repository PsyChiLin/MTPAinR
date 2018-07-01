######################### Permutation Test : Maximum Statistic (pt_ms) #########################

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

### Start MUA_pt_ms

### Generate null distribution 
permdta <- dta
# Resapmpling
permnum <- 10000
# Create a matrix to store Maximum Statistic: Maximum t value (abs max t)
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
  # Store the Maximum Statistic: Maximum t value into the matrix 
  permrst[j,1] <- max(abs(rdrst[1,]), na.rm = T)
}

###  Save the result to results folder as "MUA_pt_ms_Rst.Rdata"
# saveRDS(permrst,"Results/MUA_pt_ms_Rst.Rdata") # save the results
MUA_pt_ms_Rst <- readRDS("Results/MUA_pt_ms_Rst.Rdata")

### Show the Maximum Statistic and the significance threshold
maxstat_t <- quantile(MUA_pt_ms_Rst[,1], probs = 0.95) # or 0.975 in two tails
maxstat_t

### Results: Real data versus Significance threshold
# Read MUA results
MUA_Rst <- readRDS("Results/MUA_Rst.Rdata")

# Is there any significant time point ?
MUA_Rst$tvalue[abs(MUA_Rst$tvalue) > maxstat_t]

# The biggest t value in the real dataset is 3.025338 (smaller than maxstat_t, 3.72)  
# max(abs(MUA_Rst$tvalue)) 

### Plot the results quickly
## Note : the absolute t value of real data is used.
ggplot(data = MUA_Rst, aes(x = Times, y = abs(tvalue)))+
    geom_line(size = 1.2)+
    geom_hline(yintercept = maxstat_t, col = "red", size = 1)+
    theme_bw()+
    ylab("t value")+
    xlab("Time(s)")
