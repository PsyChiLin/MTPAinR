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
# Re-define the label of the condition
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
MUA_Rst <- melt(MUA_Rst[,c(7,8,9,10)],id.vars = "Times")
colnames(MUA_Rst) <- c("Times","Method","Adj_p")
MUA_Rst$Method <- factor(MUA_Rst$Method, levels = c("BONpv_log","BHpv_log","BYpv_log"))
Figure3_1 <- ggplot(data = MUA_Rst, aes(x = Times, y = Adj_p, col = Method))+
  geom_line(size = 1.2)+
  facet_grid(~Method)+
  geom_hline(yintercept = (-2*log10(0.05)), col = "red",size = 1)+
  theme_bw()+
  ylab("-2log10(pvalue)")+
  xlab("Time(s)")+
  scale_colour_manual(labels = c("Bonferroni","FDR(BH Method)","FDR(BY Method)"),
                      values = c("chartreuse4", "firebrick", "gold"))+
  ggtitle("(A) p-value Corrections")+
  ylim(0, 3.5)+
  theme(#legend.position = c(0.9,0.87),
        strip.background  = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0))

MUA_pt_ms_Rst <- as.data.frame(readRDS("Results/MUA_pt_ms_Rst.Rdata"))
colnames(MUA_pt_ms_Rst) <- "MUA_pt_ms"
maxstat_t <- quantile(MUA_pt_ms_Rst[,1], probs = 0.95)
MUA_pt_1dtc_Rst <- as.data.frame(readRDS("Results/MUA_pt_1dtc_Rst.Rdata"))
colnames(MUA_pt_1dtc_Rst) <- "MUA_pt_1dtc"
maxstat_stcz <- quantile(MUA_pt_1dtc_Rst[,1], probs = 0.95)

Figure3_2_1 <- ggplot(data = MUA_pt_ms_Rst, aes(x = MUA_pt_ms))+  
  geom_histogram(alpha=.5, col = "#000000", fill = "#000000" )+
  geom_vline(aes(xintercept=maxstat_t),linetype="dashed",size=1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 10))+
  ylab("Count")+
  xlab("Maximum t")+
  ggtitle("Maximum t-statistic")

Figure3_3 <- ggplot(data = MUA_pt_1dtc_Rst , aes(x = MUA_pt_1dtc))+  
  geom_histogram(alpha=.7, col = "#000000", fill = "#000000" )+
  geom_vline(aes(xintercept=maxstat_stcz),linetype="dashed",size=1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 10))+
  ylab("Count")+
  xlab("Maximum STCZ")+
  ggtitle("Maximum STCZ-statistic")

# MUA_Rst <- readRDS("Results/MUA_Rst.Rdata")
# Figure3_2_2 <- ggplot(data = MUA_Rst, aes(x = Times, y = abs(tvalue)))+
#   geom_line(size = 1.2)+
#   geom_hline(yintercept = maxstat_t,col = "red", size = 1)+
#   theme_bw()+
#   ylab("Absolute t value")+
#   xlab("Time(s)")+
#   ylim(0, 3.5)+
#   theme(plot.title = element_text(hjust = 0.5))
# Figure3_2 <- grid.arrange(Figure3_2_1,Figure3_2_2,ncol = 2, 
#                           top = textGrob("(B) Maximum t-statistic",
#                                          hjust = 2))

Figure3_23 <- grid.arrange(Figure3_2_1,Figure3_3,ncol = 2,widths = c(5,5),
                           top = textGrob("(B) Non-Parametric Permutation Frameworks",hjust = 1.7))
Figure3 <- grid.arrange(Figure3_1,Figure3_23,ncol = 1, top = "Mass Univariate Analysis")

pdf(file = "Figures/Figure3.pdf", height = 8, width = 12)
grid.draw(Figure3)
dev.off()

### Figure 4
MTPA_Rst <- readRDS("Results/MTPA_bin2_Rst.Rdata")
MAU_Rst<- readRDS("Results/MAU_GLM_Rst.Rdata")
MTPA_Rst <-  filter(MTPA_Rst, Area == "LIFG")
MTPA_Rst$Method <- "Multi Time Points Analysis"
MAU_Rst <-  filter(MAU_Rst,Area == "LIFG")
MAU_Rst$Method <- "Mass Univariate Analysis"
Rst <- rbind(MTPA_Rst,MAU_Rst)
Figure4 <- ggplot(data = Rst ,aes(x =Times, y = AUC))+
  geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Method),
              alpha = 0.3)+
  geom_line(aes(col = Method),size = 1.2)+
  theme_bw()+
  facet_grid(~Method)+
  ylab("AUC")+
  xlab("Time(s)")+
  geom_hline(yintercept = 0.5, col = "red", size = 1)+
  scale_colour_manual(values=c("chartreuse4", "firebrick"))+
  scale_fill_manual(values=c("chartreuse4", "firebrick"))+
  theme(legend.position = "none",
        strip.background  = element_blank(),
        plot.title = element_text(hjust = 0.2))

pdf(file = "Figures/Figure4.pdf", height = 5, width = 7)
Figure4
dev.off()

### Figure 5
MTPA_Rst <- readRDS("Results/MTPA_bin2_Rst.Rdata")
Figure5_1 <- ggplot(data = MTPA_Rst, aes(x =Times, y = AUC))+
  geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Area),alpha = 0.3)+
  geom_line(size = 1.2, aes(col = Area))+
  theme_bw()+
  theme(strip.background = element_rect(color = "white", fill = "white", size = 0.2),
        plot.title = element_text(hjust = 0, size = 10), legend.position = "none")+
  facet_grid(~Area)+
  ylab("AUC")+
  xlab("Time(s)")+
  ggtitle("(A) AUC with Confidence Intervals")+
  geom_hline(yintercept = 0.5, col = "red", size = 1)+
  scale_colour_manual(values=c("chartreuse4", "firebrick"))+
  scale_fill_manual(values=c("chartreuse4", "firebrick"))

Figure5_2 <- ggplot(data = MTPA_Rst, aes(x =Times, y = AUC, fill = Area))+
  geom_line(size = 1.2, aes(col = Area))+
  theme_bw()+
  ggtitle("(B) Average AUC Curves")+
  theme(plot.title = element_text(hjust = 0,size = 10))+
  ylab("AUC")+
  xlab("Time(s)")+
  scale_colour_manual(values=c("chartreuse4", "firebrick"))+
  scale_fill_manual(values=c("chartreuse4", "firebrick"))

Figure5 <- grid.arrange(Figure5_1, Figure5_2, ncol = 1, top = "Areas Comparison")

pdf(file = "Figures/Figure5.pdf", height = 7, width = 7)
grid.draw(Figure5)
dev.off()

### Supplementary Figure 1
b2 <- readRDS("Results/MTPA_bin2_Rst.Rdata")
b2 <-  filter(b2,Area == "LIFG")
b2$Bandwidth <- "Bandwidth = 2"
b3 <- readRDS("Results/MTPA_bin3_Rst.Rdata")
b3 <-  filter(b3,Area == "LIFG")   
b3$Bandwidth <- "Bandwidth = 3"
Bandwith2vs3 <- rbind(b2,b3)
SupplementaryFigure1 <- ggplot(data = Bandwith2vs3 ,aes(x =Times, y = AUC))+
  geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u, fill = Bandwidth),alpha = 0.3)+
  geom_line(aes(col = Bandwidth),size = 1.2)+
  theme_bw()+
  facet_grid(~Bandwidth)+
  ylab("AUC")+
  xlab("Time(s)")+
  geom_hline(yintercept = 0.5, col = "red", size = 1)+
  scale_colour_manual(values=c("chartreuse4", "firebrick"))+
  scale_fill_manual(values=c("chartreuse4", "firebrick"))+
  theme(legend.position = "none",strip.background  = element_blank(),plot.title = element_text(hjust = 0.2))

pdf(file = "Figures/SupplementaryFigure1.pdf", height = 4, width = 5)
SupplementaryFigure1
dev.off()

### Github README.Rmd Figure 1
MUA_Rst <- readRDS("Results/MUA_Rst.Rdata")
MUA_Rst <- melt(MUA_Rst[,c(7,8,10)],id.vars = "Times")
colnames(MUA_Rst) <- c("Times","Method","Adj_p")
MUA_Rst$Method <- factor(MUA_Rst$Method, levels = c("BONpv_log","BHpv_log"))
MTPA_Rst <- readRDS("Results/MTPA_bin2_Rst.Rdata")

GithubFigure1_1 <- ggplot(data = MUA_Rst, aes(x = Times, y = Adj_p, col = Method))+
  geom_line(size = 1.2)+
  geom_hline(yintercept = (-2*log10(0.05)), col = "red",size = 1)+
  theme_bw()+
  ylab("-2log10(pvalue)")+
  xlab("Time(s)")+
  scale_colour_manual(labels = c("Bonferroni","FDR(BH Method)"),
                      values=c("chartreuse4", "firebrick"))+
  ggtitle("MUA")+
  ylim(0, 3.5)+
  theme(legend.position = c(0.7,0.87),plot.title = element_text(hjust = 0.5, size = 12))
GithubFigure1_2 <- ggplot(data = filter(MTPA_Rst, Area == "LIFG"), aes(x =Times, y = AUC))+
  geom_ribbon(aes(ymax = AUC_l,ymin = AUC_u),fill = "firebrick",alpha = 0.3)+
  geom_line(size = 1.2, col = "firebrick")+
  theme_bw()+
  theme(strip.background = element_rect(color = "white", fill = "white", size = 0.2),
        plot.title = element_text(hjust = 0.5, size = 12), legend.position = "none")+
  ylab("AUC")+
  xlab("Time(s)")+
  ggtitle("MTPA")+
  geom_hline(yintercept = 0.5, col = "red", size = 1)


GithubFigure1 <- grid.arrange(GithubFigure1_1, GithubFigure1_2, ncol = 2)

pdf(file = "Figures/GithubFigure1.pdf", height = 5, width = 7)
grid.draw(GithubFigure1)
dev.off()