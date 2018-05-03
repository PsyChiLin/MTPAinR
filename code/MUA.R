######################### Mass Univariate Analysis (MUA) #########################

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

### Start MUA
# Define a matrix to store the results from all time points
muarst <- matrix(NA,5,167)
# Using the for loop to examine the difference at each sampled time point
for (i in 1:167){
        # the present demonstartion 
        testdta <- filter(dta[,c(1,2,3,i+3)]) 
        # fit a paired T model
        test <- t.test(testdta[,4]~testdta[,3],paired = T) 
        # Store the t value and p value into the matrix
        muarst[1,i] <- test$statistic
        muarst[2,i] <- test$p.value
}

### P value adjustment
# Bonferroni corrections
muarst[3,] <- p.adjust(muarst[2,],method = "bonferroni")
# Benjamini & Hochberg (1995) FDR corrections
muarst[4,] <- p.adjust(muarst[2,],method = "BH")
# Benjamini & Yekutieli (2001) BY-FDR corrections
muarst[5,] <- p.adjust(muarst[2,],method = "BY")

### Transform the data to a dataframe
muarst <- as.data.frame(muarst)
colnames(muarst) <- colnames(dta)[4:170]
row.names(muarst) <- c("tvalue","pvalue","BONpv","BHpv","BYpv")
muarst <- as.data.frame(t(muarst))
muarst$pvalue_log <- -(2*log10(muarst$pvalue))
muarst$BHpv_log <- -(2*log10(muarst$BHpv))
muarst$BONpv_log <- -(2*log10(muarst$BONpv))
muarst$BYpv_log <- -(2*log10(muarst$BYpv))
muarst$Times <- tp

### Save the result to results folder as "MUA_Rst.Rdata"
#saveRDS(muarst,"Results/MUA_Rst.Rdata")

### Plot the results quickly
# plt <- melt(muarst[,-c(1:5)],id.vars = "Times")
# colnames(plt)[2:3] <- c("Method","Adj_p")
# ggplot(data = plt, aes(x = Times, y = Adj_p, col = Method ))+
#         geom_line()+
#         geom_hline(yintercept = (-2*log10(0.05)), col = "red")+
#         theme_bw()+
#         ylab("-2log10(pvalue)")



