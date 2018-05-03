################# Mass Univariate Analysis (MUA) #################

## Load Supporting Tools
library(pacman)
pacman::p_load(ERP, mnormt, fdrtool, tidyverse, gridExtra, crayon,
               boot, reshape2, ggthemes, devtools,randomForest,leaps,pROC)
## Read Data
dta <- readRDS("Data/NIRSdata_LTFGLMTG.Rdata")
tp <- seq(0,16,by=0.0959)

## Explore Data
head(dta)

## Start MUA
# Define a matrix to store the results from all time points
muarst <- matrix(NA,5,178) 
# Using the for loop to 
for (i in 12:178){
        testdta <- filter(dta[,c(1,2,3,i+3)], Area == "LIFG")
        test <- t.test(testdta[,4]~testdta[,3],paired = T)
        muarst[1,i] <- test$statistic
        muarst[2,i] <- test$p.value
}
muarst[3,] <- p.adjust(muarst[2,],method = "bonferroni")
muarst[4,] <- (p.adjust(muarst[2,],method = "BH"))
muarst[5,] <- (p.adjust(muarst[2,],method = "BY"))
muarst <- as.data.frame(muarst)
colnames(muarst) <- colnames(dta)[4:181]
row.names(muarst) <- c("tvalue","pvalue","BONpv","BHpv","BYpv")

pl <- as.data.frame(t(muarst))
pl$Pvalue <- -(2*log10(pl$pvalue))
pl$BH_Pvalue <- -(2*log10(pl$BHpv))
pl$BON_Pvalue <- -(2*log10(pl$BONpv))
pl$BY_Pvalue <- -(2*log10(pl$BYpv))
pl$Times <- seq(-1,16,by=0.0958)
pl2 <- melt(pl[,-c(1:5)],id.vars = "Times")
colnames(pl2)[2:3] <- c("Method","Adj_p")

#saveRDS(pl,"mua.Rdata")
#saveRDS(pl2,"Results/mua_plot_3.Rdata")

ggplot(data = pl2, aes(x = Times, y = Adj_p, col = Method ))+
        geom_line()+
        geom_hline(yintercept = (-2*log10(0.05)), col = "red")+
        theme_bw()+
        ylab("-2log10(pvalue)")+
        ylim(0,4.67)
        #scale_colour_manual(values=c("#CC6666", "#9999CC"))



