######################### Figures #########################

### Load Supporting Tools
# clear the working directory
rm(list = ls())
# Use pacman package to manage all the packages
library(pacman)
pacman::p_load(ERP, mnormt, fdrtool,
               tidyverse, gridExtra, crayon, 
               boot, reshape2, ggthemes, grid,
               devtools,randomForest,leaps, pROC)
source("Functions/plot_tete.R")
source("Functions/g_legend.R")

### Read LIFG Data
# Read the used Data, which is collected by the present NIRS experiment, in the present study
dta <- readRDS("Data/NIRSdata_LTFGLMTG.Rdata")
# Define time point
tp <- seq(0,16,by=0.0959)
# Re-define the Area of the condition
dta$Condition <- as.factor(ifelse(dta$Condition == "0","UnRelated","Related"))

### Figure 2
Figure2_1 <- plot_tete(data = filter(dta,Area == "LIFG"),
                       channel = 1,
                       subject = 2,
                       uV = 4:170,
                       test = 3,
                       frames = tp,
                       mode = "raw",
                       labs = list(x = "Time(s)", y = "dConcentratio(uMolar)", 
                                   title = "(A) Raw Curves"),
                       ylim = c(-1e-05, 7e-06))
Figure2_2 <- plot_tete(data =  filter(dta,Area == "LIFG"),
                   channel = 1,
                   subject = 2,
                   uV = 4:170,
                   frames = tp,
                   test = 3,
                   mode = "bootci",
                   labs = list(x = "Time(s)", y = "dConcentratio(Molar)", 
                               title = "(B)  Confidence Intervals"),
                   ylim = c(-1e-05, 7e-06))
Figure2 <- grid.arrange(Figure2_1+theme(strip.text = element_blank(),
                                        legend.position = "none",
                                        axis.title.y  = element_blank(),
                                        axis.title.x  = element_blank(),
                                        plot.title = element_text(hjust = 0.1)),
                        Figure2_2+theme(strip.text = element_blank(),
                                        legend.position = c(0.7,0.2),
                                        axis.text.y = element_blank(),
                                        axis.title.y  = element_blank(),
                                        axis.title.x  = element_blank(),
                                        plot.title = element_text(hjust = 0.1)),
                        ncol = 2,
                        widths = c(2.3,2),
                        left = textGrob("dConcentratio(Molar)",rot = 90, gp=gpar(fontsize=12,font=1)),
                        bottom = textGrob("Time(s)",gp=gpar(fontsize=12,font=1)))

pdf(file = "Figures/Figure2.pdf", height = 5, width = 7)
grid.draw(Figure2)
dev.off()

### Figure 3
MUA_Rst <- readRDS("Results/MUA_Rst.Rdata")
MUA_Rst <- melt(MUA_Rst[,c(7,8,10)],id.vars = "Times")
colnames(MUA_Rst) <- c("Times","Method","Adj_p")
MUA_Rst$Method <- factor(MUA_Rst$Method, levels = c("BONpv_log","BHpv_log"))
Figure3_1 <- ggplot(data = MUA_Rst, aes(x = Times, y = Adj_p, col = Method))+
  geom_line(size = 1.2)+
  geom_hline(yintercept = (-2*log10(0.05)), col = "red",size = 1)+
  theme_bw()+
  ylab("-2log10(pvalue)")+
  xlab("Time(s)")+
  scale_colour_manual(labels = c("Bonferroni","FDR(BH Method)"),
                      values=c("chartreuse4", "firebrick"))+
  ggtitle("(A) p-value Corrections")+
  ylim(0, 3.5)+
  theme(legend.position = c(0.7,0.87),plot.title = element_text(hjust = 0.1))

permrst <- readRDS("Results/MUA_pt_ms_rst.Rdata")
maxstat_t <- quantile(permrst[,1], probs = 0.95)
MUA_Rst <- readRDS("Results/MUA_Rst.Rdata")
Figure3_2 <- ggplot(data = MUA_Rst, aes(x = Times, y = abs(tvalue)))+
  geom_line(size = 1.2)+
  geom_hline(yintercept = maxstat_t, col = "red", size = 1)+
  theme_bw()+
  ylab("Absolute t value")+
  xlab("Time(s)")+
  ggtitle("(B) Maximum Statistics")+
  ylim(0, 3.5)+
  theme(plot.title = element_text(hjust = 0.1))

Figure3 <- grid.arrange(Figure3_1,Figure3_2,ncol = 2)

pdf(file = "Figures/Figure3.pdf", height = 5, width = 7)
grid.draw(Figure3)
dev.off()


# dta3 <- readRDS("Results/AREA2_v2_Rcv100_ci90_Subterm_Default.Rdata")
# plotdta <- readRDS("Results/Area2_v2_RCV100_bin1_GLM.Rdata")
# dd1 <-  filter(dta3,Area == "LIFG")
# #dim(filter(dd1, AUC_l >0.5))
# #subset(dd1,dd1$AUC_l >0.5)
# #dim(filter(filter(dta3,Area == "LMTG"), AUC_l >0.5))
# dd1$Method <- "Multi Time Points Analysis"
# dd2 <-  filter(plotdta,Area == "LIFG")   
# dd2$Method <- "Mass Univariate Analysis"
# dd <- rbind(dd1,dd2)
# Fig5 <- ggplot(data = dd ,
#                aes(x =Times, y = AUC))+
#         geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Method),
#                     alpha = 0.3)+
#         geom_line(aes(col = Method),size = 1.2)+
#         theme_bw()+
#         facet_grid(~Method)+
#         ylab("AUC")+
#         xlab("Time(s)")+
#         geom_hline(yintercept = 0.5, col = "red", size = 1)+
#         #geom_hline(yintercept = 0.6, col = "red", size = 1, alpha = 0.5)+
#         #geom_hline(yintercept = 0.7, col = "red", size = 1, alpha = 0.7)+
#         #ggtitle("(A) RandomForest Method")+
#         scale_colour_manual(values=c("chartreuse4", "firebrick"))+
#         scale_fill_manual(values=c("chartreuse4", "firebrick"))+
#         theme(legend.position = "none",
#               strip.background  = element_blank(),
#               plot.title = element_text(hjust = 0.2))
# 
# 
# #Fig5 <- ggplot(data = filter(dta3,Area == "LIFG"),
# #       aes(x =Times, y = AUC))+
# #        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u),
# #                    alpha = 0.3)+
# #        geom_line(size = 1.2)+
# #        theme_bw()+
# #        ylab("AUC")+
# #        xlab("Time(s)")+
# #        geom_hline(yintercept = 0.5, col = "red", size = 1)+
# #        #geom_hline(yintercept = 0.6, col = "red", size = 1, alpha = 0.5)+
# #        #geom_hline(yintercept = 0.7, col = "red", size = 1, alpha = 0.7)+
# #        ggtitle("(A) RandomForest Method")+
# #        theme(plot.title = element_text(hjust = 0.2))
# 
# 
# #S1 <- ggplot(data = filter(plotdta, Area == "LIFG"),aes(x =Times, y = AUC))+
# #        geom_line(size = 1.2)+
# #        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u),
# #                    alpha = 0.3)+
# #        theme_bw()+
# #        ggtitle("(B) Mass Univariate Analysis")+
# #        theme(plot.title = element_text(hjust = 0.2),
# #              legend.position = "none")+
# #        geom_hline(yintercept = 0.5, col = "red", size = 1)+
# #        ylab("AUC")+
# #        xlab("Time(s)")
# #Fig02 <- grid.arrange(Fig5,S1, ncol =2)
# #plot(Fig02)
# 
# 
# ### Two areas Comparison
# dta3 <- readRDS("Results/AREA2_v2_Rcv100_ci90_Subterm_Default.Rdata")
# 
# mtg <- filter(dta3,Area == "LMTG") 
# mtgt1 <-  filter(dta3, AUC_l > 0.5, Area == "LMTG") # 7.9094
# ac <- filter(dta3, AUC_l > 0.5, Times < 7.9094, Area == "LIFG")
# # 27 time points IFG faster than MTG
# 
# #Fig6 <- ggplot(data =dta3,aes(x =Times, y = AUC))+
# #        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area),
# #                    alpha = 0.3)+
# #        geom_line(size = 1.2, aes(col = Area))+
# #        theme_bw()+
# #        theme(strip.background = element_rect(color = "white", fill = "white", size = 0.2),
# #              plot.title = element_text(hjust = 0.5),
# #              legend.position = "right")+
# #        #facet_grid(~Area)+
# #        ylab("AUC")+
# #        xlab("Time(s)")+
# #        ggtitle("Areas Comparison")+
# #        geom_hline(yintercept = 0.5, col = "red", size = 1)+
# #        #geom_hline(yintercept = 0.6, col = "gold", size = 1, alpha = 0.5)+
# #        #geom_hline(yintercept = 0.7, col = "gold", size = 1, alpha = 0.7)+
# #        scale_colour_manual(values=c("chartreuse4", "firebrick"))+
# #       scale_fill_manual(values=c("chartreuse4", "firebrick"))
# #        #scale_colour_manual(values=c("chartreuse4", "firebrick"))+
# #        #scale_fill_manual(values=c("chartreuse4", "firebrick"))
# #Fig03 <- Fig6 
# 
# Fig7 <- ggplot(data =dta3,aes(x =Times, y = AUC))+
#         geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area),
#                     alpha = 0.3)+
#         geom_line(size = 1.2, aes(col = Area))+
#         theme_bw()+
#         theme(strip.background = element_rect(color = "white", fill = "white", size = 0.2),
#               plot.title = element_text(hjust = 0, size = 10),
#               legend.position = "none")+
#         facet_grid(~Area)+
#         ylab("AUC")+
#         xlab("Time(s)")+
#         ggtitle("(A) AUC with Confidence Intervals")+
#         geom_hline(yintercept = 0.5, col = "red", size = 1)+
#         #geom_hline(yintercept = 0.6, col = "gold", size = 1, alpha = 0.5)+
#         #geom_hline(yintercept = 0.7, col = "gold", size = 1, alpha = 0.7)+
#         scale_colour_manual(values=c("chartreuse4", "firebrick"))+
#         scale_fill_manual(values=c("chartreuse4", "firebrick"))
# 
# #scale_colour_manual(values=c("chartreuse4", "firebrick"))+
# #scale_fill_manual(values=c("chartreuse4", "firebrick"))
# Fig8 <- ggplot(data =dta3,aes(x =Times, y = AUC, fill = Area))+
#         #geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area),
#         #            alpha = 0.3)+
#         geom_line(size = 1.2, aes(col = Area))+
#         theme_bw()+
#         ggtitle("(B) Average AUC Curves")+
#         theme(plot.title = element_text(hjust = 0,size = 10))+
#         #facet_grid(~Area)+
#         ylab("AUC")+
#         xlab("Time(s)")+
#         scale_colour_manual(values=c("chartreuse4", "firebrick"))+
#         scale_fill_manual(values=c("chartreuse4", "firebrick"))
# Fig03_2 <- grid.arrange(Fig7,Fig8,ncol = 1, 
#                       top = "Areas Comparison")
# 
# 
# 
# ### Supp 1 & 2
# #plotdta <- readRDS("Results/Area2_RCV1000_bin1_GLM.Rdata")
# #S1 <- ggplot(data = filter(plotdta, Area == "LIFG"),aes(x =Times, y = AUC))+
# #        geom_line(size = 1.2)+
# #        geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u),
# #                    alpha = 0.3)+
# #        theme_bw()+
# #        ggtitle(" Mass Univariate Analysis")+
# #        theme(plot.title = element_text(hjust = 0.5),
# #              legend.position = "none")+
# #        geom_hline(yintercept = 0.5, col = "red", size = 1)+
# #        ylab("AUC")+
# #        xlab("Time(s)")


### Supplementary Figure 1
MUA_Rst <- readRDS("Results/MUA_Rst.Rdata")
MUA_Rst <- melt(MUA_Rst[,c(6:10)],id.vars = "Times")
colnames(MUA_Rst) <- c("Times","Method","Adj_p")
MUA_Rst <- filter(MUA_Rst, Method == "BYpv_log")
SupplementaryFigure1 <- ggplot(data = MUA_Rst , aes(x = Times, y = Adj_p))+
  geom_line(size = 1.2)+
  geom_hline(yintercept = (-2*log10(0.05)), col = "red",size = 1)+
  theme_bw()+
  ylab("-2log10(pvalue)")+
  xlab("Time(s)")+
  ylim(0,3.5)+
  ggtitle("Mass Univariate Analysis: FDR(BY Method)")+
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = "Figures/SupplementaryFigure1.pdf", height = 5, width = 5)
SupplementaryFigure1
dev.off()


# b2 <- readRDS("Results/AREA2_v2_Rcv100_ci90_Subterm_Default.Rdata")
# b3 <- readRDS("Results/AREA2_v2_Rcv100_bin3_ci90_Subterm_Default.Rdata")
# b2 <-  filter(b2,Area == "LIFG")
# #dim(filter(dd1, AUC_l >0.5))
# #subset(dd1,dd1$AUC_l >0.5)
# #dim(filter(filter(dta3,Area == "LMTG"), AUC_l >0.5))
# b2$Bandwidth <- "Bandwidth = 2"
# b3 <-  filter(b3,Area == "LIFG")   
# b3$Bandwidth <- "Bandwidth = 3"
# 
# dd <- rbind(b2,b3)
# S2 <- ggplot(data = dd ,
#                aes(x =Times, y = AUC))+
#         geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Bandwidth),
#                     alpha = 0.3)+
#         geom_line(aes(col = Bandwidth),size = 1.2)+
#         theme_bw()+
#         facet_grid(~Bandwidth)+
#         ylab("AUC")+
#         xlab("Time(s)")+
#         geom_hline(yintercept = 0.5, col = "red", size = 1)+
#         #geom_hline(yintercept = 0.6, col = "red", size = 1, alpha = 0.5)+
#         #geom_hline(yintercept = 0.7, col = "red", size = 1, alpha = 0.7)+
#         #ggtitle("(A) RandomForest Method")+
#         scale_colour_manual(values=c("chartreuse4", "firebrick"))+
#         scale_fill_manual(values=c("chartreuse4", "firebrick"))+
#         theme(legend.position = "none",
#               strip.background  = element_blank(),
#               plot.title = element_text(hjust = 0.2))
# 
# 
