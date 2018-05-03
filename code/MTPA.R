library(pacman)
pacman::p_load(ERP, mnormt, fdrtool, tidyverse, gridExtra, crayon,
               boot, reshape2, ggthemes, devtools,randomForest,leaps,pROC)
setwd("~/Dropbox/R_wd/NIRS_task_ML")
dta <- readRDS("Data/AREA2_v2.Rdata")
dta <- na.omit(dta)
#head(dta)
#dim(dta)
#str(dta)
LIFG <- filter(dta, label == "LIFG")
LMTG <- filter(dta, label == "LMTG")
rm(dta)


### Set parameters
binwidth = 3 # if change binwidth, some code need to be modified
###
#nt = 10001
#ns = 7
#mt = 1
rcvnum <- 100
ci <- c(0.05,0.95)
###
upperbound <- 178-binwidth+1
lowerbound <- 12-binwidth+1
rst_LIFG <- matrix(NA,6,upperbound)
rst_LMTG <- matrix(NA,6,upperbound)

# LIFG
for (i in lowerbound:upperbound){ # change with binwidth
        if (i %in% c(seq(10,170,by =10))) {print(i)}
        ceauc <- matrix(NA,4,rcvnum)
        for (k in 1:rcvnum){
                set.seed(k)
                idc_test <- c(sample(1:14,5),sample(15:28,5))
                idc_train <- -idc_test
                fit <- randomForest(Condition~.,
                                    data = LIFG[idc_train,c(2,3,(i+3):(i+3+binwidth-1))],
                                    importance = F)
                #yhat_train_prob <- predict(fit,newdata = LIFG[idc_train,],type = "prob")[,2]
                #yhat_train_class <- predict(fit,newdata = LIFG[idc_train,],type = "class")
                yhat_test_prob <- predict(fit,newdata = LIFG[idc_test,],type = "prob")[,2]
                yhat_test_class <- predict(fit,newdata = LIFG[idc_test,],type = "class")
                #ce_train <- mean(yhat_train_class!=LIFG[idc_train,]$Condition)
                ce_test <- mean(yhat_test_class!=LIFG[idc_test,]$Condition)
                #auc_train <- pROC::auc(LIFG[idc_train,]$Condition, yhat_train_prob)
                auc_test <- pROC::auc(LIFG[idc_test,]$Condition,yhat_test_prob)
                #ceauc[1,k] <- ce_train
                ceauc[2,k] <- ce_test
                #ceauc[3,k] <- auc_train
                ceauc[4,k] <- auc_test
                
        }
        #mean(ceauc[1,])
        rst_LIFG[1,i] <- mean(ceauc[2,])
        #mean(ceauc[3,])
        rst_LIFG[2,i] <- mean(ceauc[4,])
        #quantile(ceauc[1,],probs = c(0.025,0.975))
        rst_LIFG[3,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[1]
        rst_LIFG[4,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[2]
        #quantile(ceauc[3,],probs = c(0.025,0.975))
        rst_LIFG[5,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[1]
        rst_LIFG[6,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[2]
}
LIFGm <- matrix(NA,6,length(12:upperbound))
for (i in 1:length(12:upperbound)){
        LIFGm[1,i] <- mean(rst_LIFG[1,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LIFGm[2,i] <- mean(rst_LIFG[2,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LIFGm[3,i] <- mean(rst_LIFG[3,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LIFGm[4,i] <- mean(rst_LIFG[4,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LIFGm[5,i] <- mean(rst_LIFG[5,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LIFGm[6,i] <- mean(rst_LIFG[6,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
}
LIFGm <- as.data.frame(LIFGm)
colnames(LIFGm) <- paste0("Time",1:length(12:upperbound))
row.names(LIFGm) <- c("CE","AUC","CE_l","CE_u","AUC_l","AUC_u")

# LMTG
for (i in lowerbound:upperbound){
        if (i %in% c(seq(10,170,by =10))) {print(i)}
        ceauc <- matrix(NA,4,rcvnum)
        for (k in 1:rcvnum){
                set.seed(k)
                idc_test <- c(sample(1:14,5),sample(15:28,5))
                idc_train <- -idc_test
                fit <- randomForest(Condition~.,
                                    data = LMTG[idc_train,c(2,3,(i+3):(i+3+binwidth-1))],
                                    importance = F)
                #yhat_train_prob <- predict(fit,newdata = LMTG[idc_train,],type = "prob")[,2]
                #yhat_train_class <- predict(fit,newdata = LMTG[idc_train,],type = "class")
                yhat_test_prob <- predict(fit,newdata = LMTG[idc_test,],type = "prob")[,2]
                yhat_test_class <- predict(fit,newdata = LMTG[idc_test,],type = "class")
                #ce_train <- mean(yhat_train_class!=LMTG[idc_train,]$Condition)
                ce_test <- mean(yhat_test_class!=LMTG[idc_test,]$Condition)
                #auc_train <- pROC::auc(LMTG[idc_train,]$Condition, yhat_train_prob)
                auc_test <- pROC::auc(LMTG[idc_test,]$Condition,yhat_test_prob)
                #ceauc[1,k] <- ce_train
                ceauc[2,k] <- ce_test
                #ceauc[3,k] <- auc_train
                ceauc[4,k] <- auc_test
        }
        #mean(ceauc[1,])
        rst_LMTG[1,i] <- mean(ceauc[2,])
        #mean(ceauc[3,])
        rst_LMTG[2,i] <- mean(ceauc[4,])
        #quantile(ceauc[1,],probs = c(0.025,0.975))
        rst_LMTG[3,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[1]
        rst_LMTG[4,i] <- quantile(ceauc[2,],probs = c(ci[1],ci[2]))[2]
        #quantile(ceauc[3,],probs = c(0.025,0.975))
        rst_LMTG[5,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[1]
        rst_LMTG[6,i] <- quantile(ceauc[4,],probs = c(ci[1],ci[2]))[2]
}
LMTGm <- matrix(NA,6,length(12:upperbound ))
for (i in 1:length(12:upperbound )){
        LMTGm[1,i] <- mean(rst_LMTG[1,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LMTGm[2,i] <- mean(rst_LMTG[2,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LMTGm[3,i] <- mean(rst_LMTG[3,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LMTGm[4,i] <- mean(rst_LMTG[4,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LMTGm[5,i] <- mean(rst_LMTG[5,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
        LMTGm[6,i] <- mean(rst_LMTG[6,c((i+lowerbound-1):(i+lowerbound-1+binwidth-1))])
}
LMTGm <- as.data.frame(LMTGm)
colnames(LMTGm) <- paste0("Time",1:length(12:upperbound))#
row.names(LMTGm) <- c("CE","AUC","CE_l","CE_u","AUC_l","AUC_u")

# Combine LIFG+ LMTG
temp <- as.data.frame(t(LIFGm))
temp$Area <- "LIFG"
temp$Times <- seq(-1,16,by=0.0958)[-c(1:11,upperbound+1:178)]
temp2 <- as.data.frame(t(LMTGm))
temp2$Area <- "LMTG"
temp2$Times <- seq(-1,16,by=0.0958)[-c(1:11,upperbound+1:178)]
plotdta <- rbind(temp,temp2)
#rm(temp,temp2,rst_LMTG,rst_LIFG,LIFG,LMTG,LIFGm,LMTGm,ceauc)

saveRDS(plotdta, "Results/AREA2_v2_Rcv100_bin3_ci90_Subterm_Default.Rdata")

ggplot(data = plotdta,aes(x =Times, y = AUC, col = Area))+
        geom_line()+
        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area, col = Area),
                    alpha = 0.3)+
        theme_bw()

ggplot(data = plotdta,aes(x =Times, y = AUC_l, col =Area))+
        geom_line()+
        #geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area, col = Area),
        #            alpha = 0.3)+
        theme_bw()

ggplot(data = filter(plotdta,Area == "LIFG"),aes(x =Times, y = AUC, col = Area))+
        geom_line()+
        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area, col = Area),
                    alpha = 0.3)+
        theme_bw()


ggplot(data = filter(plotdta,Area == "LMTG"),aes(x =Times, y = AUC, col = Area))+
        geom_line()+
        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area, col = Area),
                    alpha = 0.3)+
        theme_bw()
