################# permutation #################
library(pacman)
pacman::p_load(ERP, mnormt, fdrtool, tidyverse, gridExtra, crayon,dplyr,
               boot, reshape2, ggthemes, devtools,randomForest,leaps,pROC)
rm(list=ls())
setwd("~/Dropbox/R_wd/NIRS_task_ML")
dta <- readRDS("Data/AREA2_v2.Rdata")
tp <- seq(-1,16,by=0.0958)
#head(dta)

################# generate null distribution #################
permdta <- filter(dta, label == "LIFG")
permnum <- 10000
permrst <- matrix(NA,permnum,2) # t stat and p value

for (j in 1:permnum){
  print(j)
  set.seed(j) # for reproducible 
  murst <- matrix(NA,2,178) # record the statistics at each sampled time point
  permdta$Condition <- sample(permdta$Condition) # random shuffle "Condition"
  for (i in 12:178){
    testdta <- filter(permdta[,c(1,2,3,i+3)])
    test <- t.test(testdta[,4]~testdta[,3],paired = T)
    murst[1,i] <- test$statistic
    murst[2,i] <- test$p.value
  }
  permrst[j,1] <- max(murst[1,], na.rm = T)
  permrst[j,2] <- min(murst[2,], na.rm = T)
}

maxstat <- quantile(permrst[,1], probs = 0.95)
saveRDS(permrst,"Results/permutation_rst_10000.Rdata")
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
