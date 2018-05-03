
b2 <- readRDS("Results/AREA2_v2_Rcv100_ci90_Subterm_Default.Rdata")
b3 <- readRDS("Results/AREA2_v2_Rcv100_bin3_ci90_Subterm_Default.Rdata")

t.test(b2$AUC[1:165],b3$AUC[1:165],var.equal = T)
dt(-0.73,332)

mean(b2$AUC[1:165]);sd(b2$AUC[1:165])
mean(b3$AUC[1:165]);sd(b3$AUC[1:165])
?t.test



cor(b2$AUC[1:165],b3$AUC[1:165])


#t.test(b2$AUC[167:331],b3$AUC[166:330]) 

head(b3[,c(2,7,8)])
head(b2[1:165,c(2,7,8)])

cbind(b2[167:331,c(2,7,8)],b3[166:330,c(2,7,8)])

filter(b2[1:166,], AUC_l > 0.5)
filter(b3[1:165,], AUC_l > 0.5)


filter(b2[167:331,], AUC_l > 0.5)
filter(b3[166:330,], AUC_l > 0.5)
