######################### Permutation Test : Maximum Statistic (pt_ms) #########################

### Load Supporting Tools
# Use pacman package to manage all the packages
library(pacman)
pacman::p_load(ERP, mnormt, fdrtool,
               tidyverse, gridExtra, crayon, 
               boot, reshape2, ggthemes, 
               devtools,randomForest,leaps, pROC)

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
# Create a matrix to store Maximum Statistic: Maximum t value
permrst <- matrix(NA,permnum,1)

for (j in 1:permnum){
  # Print to see the progress
  print(j)
  # Set seed for reproduciblility
  set.seed(j) 
  # Record the statistics at each sampled time point
  rdrst <- matrix(NA,2,170)
  # Randomly shuffle "Condition"
  permdta$Condition <- sample(permdta$Condition) 
  # From 1 to 167 time points
  for (i in 4:170){
    # Each sampled time point
    testdta <- permdta
    # Fit a paired T model
    test <- t.test(testdta[,4]~testdta[,3],paired = T)
    # Store the t value and p value into the matrix
    rdrst[1,i] <- test$statistic
    rdrst[2,i] <- test$p.value
  }
  # Store the Maximum Statistic: Maximum t value into the matrix 
  permrst[j,1] <- max(rdrst[1,], na.rm = T)
}


maxstat <- quantile(permrst[,1], probs = 0.95)

#saveRDS(permrst,"Results/MUA_pt_ms_rst.Rdata")













################# compare with t value of real data  #################
dta2 <- readRDS("Results/mua.Rdata")
dta3 <- readRDS("Results/mua_plot_3.Rdata")
murst <- matrix(NA,5,178)
#i = 1
for (i in 12:178){
  testdta <- filter(dta[,c(1,2,3,i+3)], label == "LIFG")
  test <- t.test(testdta[,4]~testdta[,3],paired = T)
  murst[1,i] <- test$statistic
  murst[2,i] <- test$p.value
}
#murst[3,] <-(p.adjust(murst[2,],method = "BH"))

-2*log10(murst[3,])
#murst[2,] < 0.05
sigtp <- abs(murst[1,][12:178])[abs(murst[1,])[12:178] > maxstat]
length(sigtp)

realdta_abst <- na.omit(data.frame(tvalue = abs(murst[1,]),
                                   timepoint = tp))
ggplot(data = realdta_abst, 
       aes(x = timepoint, y = tvalue))+
    geom_line(size = 1.2)+
    geom_hline(yintercept = maxstat, col = "red", size = 1)+
    theme_bw()+
    ylab("t value")+
    xlab("Time(s)")


# murst[3,] <- p.adjust(murst[2,],method = "bonferroni")
# murst[4,] <- (p.adjust(murst[2,],method = "BH"))
# murst[5,] <- (p.adjust(murst[2,],method = "BY"))
# murst <- as.data.frame(murst)
# colnames(murst) <- colnames(dta)[5:181]
# row.names(murst) <- c("tvalue","pvalue","BONpv","BHpv","BYpv")
# 
# pl <- as.data.frame(t(murst))
# pl$Pvalue <- -(2*log10(pl$pvalue))
# pl$BH_Pvalue <- -(2*log10(pl$BHpv))
# pl$BON_Pvalue <- -(2*log10(pl$BONpv))
# pl$BY_Pvalue <- -(2*log10(pl$BYpv))
# pl$Times <- seq(-1,16,by=0.0958)
# pl2 <- melt(pl[,-c(1:5)],id.vars = "Times")
# colnames(pl2)[2:3] <- c("Method","Adj_p")
# 
# #saveRDS(pl,"mua.Rdata")
# #saveRDS(pl2,"Results/mua_plot_3.Rdata")
# 
# ggplot(data = pl2, aes(x = Times, y = Adj_p, col = Method ))+
#   geom_line()+
#   geom_hline(yintercept = (-2*log10(0.05)), col = "red")+
#   theme_bw()+
#   ylab("-2log10(pvalue)")+
#   ylim(0,4.67)
# #scale_colour_manual(values=c("#CC6666", "#9999CC"))
# 
# 
# 
