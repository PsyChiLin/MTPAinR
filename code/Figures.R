library(pacman)
pacman::p_load(ERP, mnormt, fdrtool, tidyverse, gridExtra, crayon,
               boot, reshape2, ggthemes, devtools,grid)
setwd("~/Dropbox/R_wd/NIRS_task_ML")
source("Functions/plot_tete.R")
source("Functions/g_legend.R")


dta <- readRDS("Data/AREA2_v2.Rdata")
dta$Condition <- as.factor(ifelse(dta$Condition == "0","UnRelated","Related"))
dta <- dta[,-c(4:14)]
tp <- seq(0,16,by=0.096)
Fig1 <- plot_tete(data = filter(dta,label == "LIFG"),
                   channel = 1,
                   subject = 2,
                   uV = 4:170,
                   frames = tp,
                   test = 3,
                   mode = "raw",
                  labs = list(x = "Time(s)", y = "dConcentratio(uMolar)", 
                              title = "(A) Raw Curves"),
                   ylim = c(-1e-05, 7e-06))
Fig2 <- plot_tete(data = filter(dta,label == "LIFG"),
                   channel = 1,
                   subject = 2,
                   uV = 4:170,
                   frames = tp,
                   test = 3,
                   mode = "mean",
                  labs = list(x = "Time(s)", y = "dConcentratio(uMolar)", 
                              title = "(B) Average Curves"),
                   ylim = c(-1e-05, 7e-06))
Fig3 <- plot_tete(data =  filter(dta,label == "LIFG"),
                   channel = 1,
                   subject = 2,
                   uV = 4:170,
                   frames = tp,
                   test = 3,
                   mode = "bootci",
                   labs = list(x = "Time(s)", y = "dConcentratio(Molar)", 
                               title = "(C)  Confidence Intervals"),
                   ylim = c(-1e-05, 7e-06))
Fig01 <- grid.arrange(Fig1+theme(strip.text = element_blank(),
                                 legend.position = "none",
                                 #axis.text.y  = element_blank(),
                                 axis.title.y  = element_blank(),
                                 axis.title.x  = element_blank(),
                                 plot.title = element_text(hjust = 0.2)),
                      Fig2+theme(strip.text = element_blank(),
                                 legend.position = "none",
                                 axis.text.y = element_blank(),
                                 axis.title.y =   element_blank(),
                                 axis.title.x  = element_blank(),
                                 plot.title = element_text(hjust = 0.2)),
                      Fig3+theme(strip.text = element_blank(),
                                 legend.position = c(0.7,0.2),
                                 axis.text.y = element_blank(),
                                 axis.title.y  = element_blank(),
                                 axis.title.x  = element_blank(),
                                 plot.title = element_text(hjust = 0.2)),
                      ncol = 3,
                      widths = c(2.3,2,2),
                      left =  textGrob("dConcentratio(Molar)",rot = 90,gp=gpar(fontsize=12,font=1)),
                      bottom = textGrob("Time(s)",gp=gpar(fontsize=12,font=1)))

dta2 <- readRDS("Results/mua_plot_3.Rdata")
dta2 <- filter(dta2, !Method %in% c("Pvalue","BY_Pvalue"))
dta2$Method <- droplevels(dta2$Method)
dta2$Method <- factor(dta2$Method, levels = c("BON_Pvalue","BH_Pvalue"))
Fig4 <- ggplot(data = dta2, aes(x = Times, y = Adj_p, col = Method))+
        geom_line(size = 1.2)+
        geom_hline(yintercept = (-2*log10(0.05)), col = "red",size = 1)+
        theme_bw()+
        ylab("-2log10(pvalue)")+
        xlab("Time(s)")+
        ylim(0,4.67)+
        scale_colour_manual(labels = c("Bonferroni","FDR(BH Method)"),
                            values=c("chartreuse4", "firebrick"))+
        ggtitle("Mass Univariate Analysis")+
        theme(legend.position = c(0.7,0.77),
              plot.title = element_text(hjust = 0.5))


dta3 <- readRDS("Results/AREA2_v2_Rcv100_ci90_Subterm_Default.Rdata")
plotdta <- readRDS("Results/Area2_v2_RCV100_bin1_GLM.Rdata")
dd1 <-  filter(dta3,Area == "LIFG")
#dim(filter(dd1, AUC_l >0.5))
#subset(dd1,dd1$AUC_l >0.5)
#dim(filter(filter(dta3,Area == "LMTG"), AUC_l >0.5))
dd1$Method <- "Multi Time Points Analysis"
dd2 <-  filter(plotdta,Area == "LIFG")   
dd2$Method <- "Mass Univariate Analysis"
dd <- rbind(dd1,dd2)
Fig5 <- ggplot(data = dd ,
               aes(x =Times, y = AUC))+
        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Method),
                    alpha = 0.3)+
        geom_line(aes(col = Method),size = 1.2)+
        theme_bw()+
        facet_grid(~Method)+
        ylab("AUC")+
        xlab("Time(s)")+
        geom_hline(yintercept = 0.5, col = "red", size = 1)+
        #geom_hline(yintercept = 0.6, col = "red", size = 1, alpha = 0.5)+
        #geom_hline(yintercept = 0.7, col = "red", size = 1, alpha = 0.7)+
        #ggtitle("(A) RandomForest Method")+
        scale_colour_manual(values=c("chartreuse4", "firebrick"))+
        scale_fill_manual(values=c("chartreuse4", "firebrick"))+
        theme(legend.position = "none",
              strip.background  = element_blank(),
              plot.title = element_text(hjust = 0.2))


#Fig5 <- ggplot(data = filter(dta3,Area == "LIFG"),
#       aes(x =Times, y = AUC))+
#        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u),
#                    alpha = 0.3)+
#        geom_line(size = 1.2)+
#        theme_bw()+
#        ylab("AUC")+
#        xlab("Time(s)")+
#        geom_hline(yintercept = 0.5, col = "red", size = 1)+
#        #geom_hline(yintercept = 0.6, col = "red", size = 1, alpha = 0.5)+
#        #geom_hline(yintercept = 0.7, col = "red", size = 1, alpha = 0.7)+
#        ggtitle("(A) RandomForest Method")+
#        theme(plot.title = element_text(hjust = 0.2))


#S1 <- ggplot(data = filter(plotdta, Area == "LIFG"),aes(x =Times, y = AUC))+
#        geom_line(size = 1.2)+
#        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u),
#                    alpha = 0.3)+
#        theme_bw()+
#        ggtitle("(B) Mass Univariate Analysis")+
#        theme(plot.title = element_text(hjust = 0.2),
#              legend.position = "none")+
#        geom_hline(yintercept = 0.5, col = "red", size = 1)+
#        ylab("AUC")+
#        xlab("Time(s)")
#Fig02 <- grid.arrange(Fig5,S1, ncol =2)
#plot(Fig02)


### Two areas Comparison
dta3 <- readRDS("Results/AREA2_v2_Rcv100_ci90_Subterm_Default.Rdata")

mtg <- filter(dta3,Area == "LMTG") 
mtgt1 <-  filter(dta3, AUC_l > 0.5, Area == "LMTG") # 7.9094
ac <- filter(dta3, AUC_l > 0.5, Times < 7.9094, Area == "LIFG")
# 27 time points IFG faster than MTG

#Fig6 <- ggplot(data =dta3,aes(x =Times, y = AUC))+
#        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area),
#                    alpha = 0.3)+
#        geom_line(size = 1.2, aes(col = Area))+
#        theme_bw()+
#        theme(strip.background = element_rect(color = "white", fill = "white", size = 0.2),
#              plot.title = element_text(hjust = 0.5),
#              legend.position = "right")+
#        #facet_grid(~Area)+
#        ylab("AUC")+
#        xlab("Time(s)")+
#        ggtitle("Areas Comparison")+
#        geom_hline(yintercept = 0.5, col = "red", size = 1)+
#        #geom_hline(yintercept = 0.6, col = "gold", size = 1, alpha = 0.5)+
#        #geom_hline(yintercept = 0.7, col = "gold", size = 1, alpha = 0.7)+
#        scale_colour_manual(values=c("chartreuse4", "firebrick"))+
#       scale_fill_manual(values=c("chartreuse4", "firebrick"))
#        #scale_colour_manual(values=c("chartreuse4", "firebrick"))+
#        #scale_fill_manual(values=c("chartreuse4", "firebrick"))
#Fig03 <- Fig6 

Fig7 <- ggplot(data =dta3,aes(x =Times, y = AUC))+
        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area),
                    alpha = 0.3)+
        geom_line(size = 1.2, aes(col = Area))+
        theme_bw()+
        theme(strip.background = element_rect(color = "white", fill = "white", size = 0.2),
              plot.title = element_text(hjust = 0, size = 10),
              legend.position = "none")+
        facet_grid(~Area)+
        ylab("AUC")+
        xlab("Time(s)")+
        ggtitle("(A) AUC with Confidence Intervals")+
        geom_hline(yintercept = 0.5, col = "red", size = 1)+
        #geom_hline(yintercept = 0.6, col = "gold", size = 1, alpha = 0.5)+
        #geom_hline(yintercept = 0.7, col = "gold", size = 1, alpha = 0.7)+
        scale_colour_manual(values=c("chartreuse4", "firebrick"))+
        scale_fill_manual(values=c("chartreuse4", "firebrick"))

#scale_colour_manual(values=c("chartreuse4", "firebrick"))+
#scale_fill_manual(values=c("chartreuse4", "firebrick"))
Fig8 <- ggplot(data =dta3,aes(x =Times, y = AUC, fill = Area))+
        #geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area),
        #            alpha = 0.3)+
        geom_line(size = 1.2, aes(col = Area))+
        theme_bw()+
        ggtitle("(B) Average AUC Curves")+
        theme(plot.title = element_text(hjust = 0,size = 10))+
        #facet_grid(~Area)+
        ylab("AUC")+
        xlab("Time(s)")+
        scale_colour_manual(values=c("chartreuse4", "firebrick"))+
        scale_fill_manual(values=c("chartreuse4", "firebrick"))
Fig03_2 <- grid.arrange(Fig7,Fig8,ncol = 1, 
                      top = "Areas Comparison")



### Supp 1 & 2
#plotdta <- readRDS("Results/Area2_RCV1000_bin1_GLM.Rdata")
#S1 <- ggplot(data = filter(plotdta, Area == "LIFG"),aes(x =Times, y = AUC))+
#        geom_line(size = 1.2)+
#        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u),
#                    alpha = 0.3)+
#        theme_bw()+
#        ggtitle(" Mass Univariate Analysis")+
#        theme(plot.title = element_text(hjust = 0.5),
#              legend.position = "none")+
#        geom_hline(yintercept = 0.5, col = "red", size = 1)+
#        ylab("AUC")+
#        xlab("Time(s)")
plotdta2 <- readRDS("Results/mua_plot_3.Rdata")
plotdta2 <- filter(plotdta2, Method == "BY_Pvalue")
#plotdta2$Method <- factor(plotdta2$Method, levels = c("BON_Pvalue","BH_Pvalue","BY_Pvalue"))
S1 <- ggplot(data = plotdta2 , aes(x = Times, y = Adj_p))+
        geom_line(size = 1.2)+
        geom_hline(yintercept = (-2*log10(0.05)), col = "red",size = 1)+
        theme_bw()+
        ylab("-2log10(pvalue)")+
        xlab("Time(s)")+
        ylim(0,4.67)+
        #scale_colour_manual(labels = c("Bonferroni","FDR(BH Method)","FDR(BY Method)"),
        #                    values=c("#CC6666", "#9999CC"))+
        ggtitle("Mass Univariate Analysis : FDR(BY Method)")+
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))

b2 <- readRDS("Results/AREA2_v2_Rcv100_ci90_Subterm_Default.Rdata")
b3 <- readRDS("Results/AREA2_v2_Rcv100_bin3_ci90_Subterm_Default.Rdata")
b2 <-  filter(b2,Area == "LIFG")
#dim(filter(dd1, AUC_l >0.5))
#subset(dd1,dd1$AUC_l >0.5)
#dim(filter(filter(dta3,Area == "LMTG"), AUC_l >0.5))
b2$Bandwidth <- "Bandwidth = 2"
b3 <-  filter(b3,Area == "LIFG")   
b3$Bandwidth <- "Bandwidth = 3"

dd <- rbind(b2,b3)
S2 <- ggplot(data = dd ,
               aes(x =Times, y = AUC))+
        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Bandwidth),
                    alpha = 0.3)+
        geom_line(aes(col = Bandwidth),size = 1.2)+
        theme_bw()+
        facet_grid(~Bandwidth)+
        ylab("AUC")+
        xlab("Time(s)")+
        geom_hline(yintercept = 0.5, col = "red", size = 1)+
        #geom_hline(yintercept = 0.6, col = "red", size = 1, alpha = 0.5)+
        #geom_hline(yintercept = 0.7, col = "red", size = 1, alpha = 0.7)+
        #ggtitle("(A) RandomForest Method")+
        scale_colour_manual(values=c("chartreuse4", "firebrick"))+
        scale_fill_manual(values=c("chartreuse4", "firebrick"))+
        theme(legend.position = "none",
              strip.background  = element_blank(),
              plot.title = element_text(hjust = 0.2))


