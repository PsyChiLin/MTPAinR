plot_tete <- function(data,  # The ERP or NIRS data
                    frames = NULL, # The time point of the data
                    uV , # The column of "Data"
                    subject = NULL, # The column of "Subject"
                    channel = NULL,  # The column of "Channel"
                    test = NULL, # The column of the variable that you want to compare (e.g. "Condition")
                    mode = c("raw","mean","bootci"), # five modes of the possible plot 
                    order = F, # ACF for scalp = FALSE ####
                    curve.col = NULL, # the color for the curve (e.g. two conditions need two colors)
                    labs = list(x = "Time (ms)", y = "Amplitude (microvolt)"),# Labs of the plot
                    ggtheme = NULL, # set your own theme
                    scalp = FALSE, # Do you want to plot on the scalp ?
                    coord.mat = NULL, # read your own coordinate matrix, default will be the full 10/10 system.
                    ylim = c(-20,20), # y limits of the plot
                    curve.fun = function(x, d){return(mean(x[d]))}, # when  mode == mean or bootci
                    boot.num = 200 ,    # when mode == bootci, the bootstraping number 
                    boot.intval = c(0.025,0.975), # when mode == bootci, the bootstraping interval
                    ci.alpha = 0.3 ){ # when mode == bootci, the alpha value of the ribbon interval
        # Some Check Functions (need More)
        options(warn=-1)
        # the function of selecting data
        dta <- data
        # Some Check Functions (need More)
        if (scalp == TRUE & order == TRUE){message("The function will not order the curve when scalp = TRUE")}
        if (is.null(frames)){frames  = 1:length(uV)}
        if (length(levels(dta[,channel])) == 1 & scalp == T){ 
                message("A single channel data could not be placed on scalp ! change to scalp = FALSE")
                scalp = FALSE
        }
        if (length(levels(dta[,channel])) == 1 & order == T){ 
                message("A single channel data could not be ordered ! change to order = FALSE ")
                order = FALSE
        }
        if (!mode %in% c("raw","mean","bootci")){stop("mode should be 'raw','mean','boot.ci'.")}
        if (length(mode) != 1) {stop("You should select only one mode of plot !")}
        ### test variable issue
        if (!is.null(test)){
                if (!class(dta[,test]) %in% c("numeric","factor")){stop("test should be numeric or factor" )}
                if (class(dta[,test]) == "factor"){
                        if (is.null(curve.col)){
                                if (length(levels(as.factor(dta[,test]))) == 1){curve.col = c("chartreuse4")}
                                if (length(levels(as.factor(dta[,test]))) == 2){curve.col = c("chartreuse4", "firebrick")}
                                if (length(levels(as.factor(dta[,test]))) == 3){curve.col = c("chartreuse4", "firebrick","gold")}
                                # add some colors
                                message("Using Default color as the curve color")
                        }
                        if (length((levels(dta[,test]))) != length(curve.col)){stop("curve.col should be equal to levels of test !")}
                }
                if (class(dta[,test]) == "numeric"){
                        if (is.null(curve.col)) {
                                if (mode == "bootci"){curve.col = c("chartreuse4")}
                                if (mode == "raw"){curve.col = c("#132B43","#56B1F7")}
                                #if (length(curve.col)!=2) {stop("curve.col should have 2 color !")}
                        }
                }
        }
        if (is.null(test)){curve.col = c("chartreuse4")}
        ## Three modes of test
        # 1. No test and test all the same
        # 2. Two or more test level
        # 3. numeric
        if (is.null(test)|length(levels(dta[,test])) == 1) {testmode = "single"}
        if (length(levels(dta[,test])) > 1) {testmode = "typical"}
        if (class(dta[,test]) == "numeric") {testmode = "numer"}
        # Some plot functions
        edaplot <- function(data,frames=NULL,uV,subject=NULL,channel=NULL){
                colnames(data)[channel] <- "Channel"
                oth.var = (1:dim(data)[2]) [! 1:dim(data)[2] %in% c(uV,subject,channel)]
                data$Info <- NA
                for (i in 1 : dim(data)[1]){
                        info1 <- paste0(as.character(data[i,oth.var]),collapse=";")
                        data$Info[i] <- paste(data[i,channel],data[i,subject],info1,sep=";")
                }
                datalong <- melt(data,
                                 id=c(variable.names(data)[c(subject,channel,oth.var)],
                                      "Info"))
                datalong <- datalong[order(datalong$Info),]
                datalong$frames <- rep(frames,length(datalong[,1])/length(frames))
                plot <- ggplot(datalong,aes(x=frames,y=value,group=Info))
                return(plot)
        }
        # The summarize function
        data_summarize <- function(data,uV,summary.var,fun=mean,...){ 
                dta <- data
                options(warn=-1)
                agglength <- length(summary.var)
                aggvar_list <- list(dta[,summary.var[1]])
                if (agglength > 1){
                        for (i in 2:agglength ){
                                aggvar_list <- append(aggvar_list,list(dta[,summary.var[2]]))
                        }
                }
                aggdata <- aggregate(dta[,uV],by=aggvar_list,
                                     fun,...)
                aggdata <- aggdata[,1:(agglength+length(uV))]
                for (i in 1: agglength){
                        colnames(aggdata)[i] <- colnames(dta)[summary.var[i]]
                }
                rownames(aggdata) <- 1:dim(aggdata)[1] 
                return(aggdata)
        }
        # The bootstraping function
        bootstrap <- function(x,bootnum,bootfun,bootintval=c(0.025,0.975),quantilenum,...){
                boot_result <- boot(x,statistic = bootfun,R = bootnum,...)
                return(quantile(boot_result$t,bootintval,na.rm = T)[quantilenum])
        }
        # Produce the data of confidence interval
        if (mode == "bootci"){
                if (testmode == "single"){
                        data_FUN <- data_summarize(dta,uV,summary.var = channel,fun = curve.fun)
                        data_Q1 <- data_summarize(dta,uV,summary.var = channel,
                                                  fun=bootstrap,bootnum=boot.num,bootfun=curve.fun,
                                                  bootintval=boot.intval,quantilenum=1)
                        data_Q2 <- data_summarize(dta,uV,summary.var = channel,
                                                  fun=bootstrap,bootnum=boot.num,bootfun=curve.fun,
                                                  bootintval=boot.intval,quantilenum=2)
                        data_fun_long <- melt(data_FUN,id=c(colnames(dta)[channel]))
                        colnames(data_fun_long)[3] <- "FUN"
                        data_Q1_long <- melt(data_Q1,id=c(colnames(dta)[channel]))
                        colnames(data_Q1_long)[3] <- "Q1"
                        data_Q2_long <- melt(data_Q2,id=c(colnames(dta)[channel]))
                        colnames(data_Q2_long)[3] <- "Q2"
                        data_ci <- full_join(data_fun_long,data_Q1_long,by = c(colnames(dta)[channel],"variable"))
                        data_ci <- full_join(data_ci,data_Q2_long,by = c(colnames(dta)[channel],"variable"))
                        colnames(data_ci)[1] <- c("Channel")
                } 
                if (testmode == "typical") {
                        dta[,test] <- as.factor(dta[,test])
                        dta[,channel] <- as.factor(dta[,channel])
                        data_FUN <- data_summarize(dta,uV,summary.var = c(channel,test),fun = curve.fun)
                        data_Q1 <- data_summarize(dta,uV,summary.var = c(channel,test),
                                                  fun=bootstrap,bootnum=boot.num,bootfun=curve.fun,
                                                  bootintval=boot.intval,quantilenum=1)
                        data_Q2 <- data_summarize(dta,uV,summary.var = c(channel,test),
                                                  fun=bootstrap,bootnum=boot.num,bootfun=curve.fun,
                                                  bootintval=boot.intval,quantilenum=2)
                        data_fun_long <- melt(data_FUN,id=c(colnames(dta)[channel],colnames(dta)[test]))
                        colnames(data_fun_long)[c(1,2,4)] <- c("Channel","Condition","FUN")
                        data_Q1_long <- melt(data_Q1,id=c(colnames(dta)[channel],colnames(dta)[test]))
                        colnames(data_Q1_long)[c(1,2,4)] <- c("Channel","Condition","Q1")
                        data_Q2_long <- melt(data_Q2,id=c(colnames(dta)[channel],colnames(dta)[test]))
                        colnames(data_Q2_long)[c(1,2,4)] <- c("Channel","Condition","Q2")
                        data_ci <- full_join(data_fun_long,data_Q1_long,by = c("Channel","Condition","variable"))
                        data_ci <- full_join(data_ci,data_Q2_long,by = c("Channel","Condition","variable"))
                } 
                #if (class(dta[,test]) == "numeric"){stop("Bootci only work when test is a factor variable.")}
        }
        # Produce default geom_theme
        if (is.null(ggtheme)){
                if (scalp == TRUE){
                        theme_default <- function(base_size = 6, base_family = ""){ 
                                theme_bw(base_size = base_size, base_family = base_family) %+replace% 
                                        theme( 
                                                panel.border     = element_blank(), 
                                                axis.line        = element_line(colour = "black"), 
                                                panel.grid.major = element_line(), 
                                                panel.grid.major.x = element_blank(), 
                                                panel.grid.major.y = element_blank(), 
                                                panel.grid.minor = element_blank(), 
                                                panel.grid.minor.x = element_blank(), 
                                                panel.grid.minor.y = element_blank(), 
                                                strip.background = element_rect(color = "white", size = 0.2), 
                                                legend.key       = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.title.x = element_blank(),
                                                axis.text.x = element_blank(),                                                                 axis.ticks =  element_blank()
                                        ) 
                        }
                }
                if (scalp == FALSE){
                        theme_default <- function(base_size = 10, base_family = ""){ 
                                theme_bw(base_size = base_size, base_family = base_family) %+replace% 
                                        theme( strip.background = element_blank()
                                        ) 
                        }
                }
        } else {theme_default <- ggtheme}
        # Draw the plots
        if (scalp == F){
                if (mode == "raw"){
                        if (order == TRUE){
                                if (testmode == "single"|testmode == "numer"){
                                        data_plot <- data_summarize(dta,uV,summary.var = channel,fun = curve.fun)
                                        colnames(data_plot)[1] <- "Channel"
                                        dta_clustorder <- data_plot
                                        hc.pred <- hclust(diss(as.matrix(dta_clustorder[,-1]), "ACF"))
                                        dta_clustorder <- dta_clustorder[hc.pred$order,]
                                        dta$Channel <- factor(dta$Channel,levels = c(as.character(dta_clustorder$Channel)))
                                }
                                if (testmode == "typical"){
                                        data_plot <- data_summarize( filter(dta,dta[,test] == levels(dta[,test])[1]),uV ,summary.var = channel,fun = curve.fun)
                                        colnames(data_plot)[1] <- "Channel"
                                        dta_clustorder <- data_plot
                                        hc.pred <- hclust(diss(as.matrix(dta_clustorder[,-1]), "ACF"))
                                        dta_clustorder <- dta_clustorder[hc.pred$order,]
                                        dta$Channel <- factor(dta$Channel,levels = c(as.character(dta_clustorder$Channel)))
                                }
                        }
                        data_plot <- dta
                        colnames(data_plot)[channel] <- "Channel"
                        if (testmode == "single"){
                                plot <- edaplot(data = data_plot,frames = frames,uV = uV,subject=subject,
                                                channel=channel)+
                                        geom_line(col = curve.col[1])+
                                        facet_wrap(~Channel)+
                                        labs(labs)+
                                        ylim(ylim)+
                                        theme_default()
                        } 
                        if (testmode == "typical"){
                                colnames(data_plot)[test] <- "Condition"
                                values <- curve.col[1]
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                plot <- edaplot(data = data_plot,frames = frames,uV = uV,subject=subject,
                                                channel=channel)+
                                        geom_line(aes(col = Condition))+
                                        facet_wrap(~Channel)+
                                        labs(labs)+
                                        scale_color_manual(values = values, name = colnames(dta)[test])+
                                        ylim(ylim)+
                                        theme_default()
                        }
                        if (testmode == "numer"){
                                colnames(data_plot)[test] <- "numer"
                                plot <- edaplot(data = data_plot,frames = frames,uV = uV,subject=subject,
                                                channel=channel)+
                                        geom_line(aes(col = numer))+
                                        facet_wrap(~Channel)+
                                        labs(labs)+
                                        scale_color_gradient(low =curve.col[1], high = curve.col[2],name =  colnames(data)[test])+
                                        ylim(ylim)+
                                        theme_default()
                         
                        }
                }
                if (mode == "mean") {
                        if (testmode == "single") {
                                data_plot <- data_summarize(dta,uV ,summary.var = channel,fun = curve.fun)
                                if (order == TRUE){
                                        colnames(data_plot)[1] <- "Channel"
                                        dta_clustorder <- data_plot
                                        hc.pred <- hclust(diss(as.matrix(dta_clustorder[,-1]), "ACF"))
                                        dta_clustorder <- dta_clustorder[hc.pred$order,]
                                        data_plot$Channel <- factor(data_plot$Channel,levels = c(as.character(dta_clustorder$Channel)))
                                }
                                data_plot <-  suppressMessages(melt(data_plot))
                                colnames(data_plot)[1] <- "Channel"
                                data_plot <- data_plot[order(data_plot$Channel,data_plot$variable),]
                                data_plot$frames <- rep(frames,dim(data_plot)[1]/length(frames))
                                plot <- ggplot(data = data_plot,aes(x = frames, y = value))+
                                        geom_line(col = curve.col[1])+
                                        facet_wrap(~Channel)+
                                        theme_bw()+
                                        labs(labs)+
                                        ylim(ylim)+
                                        theme_default()
                        } 
                        if (testmode == "typical") {
                                values <- curve.col[1]
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                data_plot <- data_summarize(dta,uV ,summary.var = c(channel,test),fun = curve.fun)
                                if (order == TRUE){
                                        colnames(data_plot)[c(1,2)] <- c("Channel","Condition")
                                        dta_clustorder <- filter(data_plot,data_plot$Condition == levels(data_plot$Condition)[1])
                                        hc.pred <- hclust(diss(as.matrix(dta_clustorder[,-c(1:2)]), "ACF"))
                                        dta_clustorder <- dta_clustorder[hc.pred$order,]
                                        data_plot$Channel <- factor(data_plot$Channel,levels = c(as.character(dta_clustorder$Channel)))
                                }
                                data_plot <-  suppressMessages(melt(data_plot))
                                colnames(data_plot)[c(1,2)] <- c("Channel","Condition")
                                data_plot <- data_plot[order(data_plot$Channel,data_plot$Condition,data_plot$variable),]
                                data_plot$frames <- rep(frames,dim(data_plot)[1]/length(frames))
                                plot <- ggplot(data = data_plot,aes(x = frames, y = value))+
                                        geom_line(aes(col = Condition))+
                                        facet_wrap(~Channel)+
                                        theme_bw()+
                                        labs(labs)+
                                        scale_color_manual(values = values, name = colnames(dta)[test])+
                                        ylim(ylim)+
                                        theme_default()
                        }
                        if (testmode == "numer") {stop("Numeric variable could not adopt mean curve ! It will be meaningless.")}
                } 
                if (mode == "bootci"){
                        if (testmode == "single"){
                                if (order == TRUE){
                                        data_plot <- data_summarize(dta,uV ,summary.var = channel,fun = curve.fun)
                                        colnames(data_plot)[1] <- "Channel"
                                        dta_clustorder <- data_plot
                                        hc.pred <- hclust(diss(as.matrix(dta_clustorder[,-1]), "ACF"))
                                        dta_clustorder <- dta_clustorder[hc.pred$order,]
                                        data_ci$Channel <- factor(data_ci$Channel,levels = c(as.character(dta_clustorder$Channel)))
                                }
                                data_ci <- data_ci[order(data_ci$Channel),]
                                data_ci$frames <- rep(frames,dim(data_ci)[1]/length(frames))
                                plot <- ggplot(data = data_ci,aes(x = frames, y = FUN))+
                                        geom_ribbon(aes(ymax = Q1,ymin = Q2),alpha = ci.alpha, fill = curve.col[1])+
                                        geom_line(col = curve.col[1])+
                                        facet_wrap(~Channel)+
                                        labs(labs)+
                                        ylim(ylim)+
                                        theme_default()
                        } 
                        if (testmode == "typical"){
                                if (order == TRUE){
                                        data_plot <- data_summarize(dta,uV ,summary.var = c(channel,test),fun = curve.fun)
                                        colnames(data_plot)[c(1,2)] <- c("Channel","Condition")
                                        dta_clustorder <- filter(data_plot,data_plot$Condition == levels(data_plot$Condition)[1])
                                        hc.pred <- hclust(diss(as.matrix(dta_clustorder[,-c(1:2)]), "ACF"))
                                        dta_clustorder <- dta_clustorder[hc.pred$order,]
                                        data_ci$Channel <- factor(data_ci$Channel,levels = c(as.character(dta_clustorder$Channel)))
                                }
                                values <- curve.col[1]
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                data_ci <- data_ci[order(data_ci$Channel,data_ci$Condition,data_ci$variable),]
                                data_ci$frames <- rep(frames,dim(data_ci)[1]/length(frames))
                                plot <- ggplot(data = data_ci,aes(x = frames, y = FUN))+
                                        geom_ribbon(aes(ymax = Q2,ymin=Q1, fill = Condition),alpha = ci.alpha)+
                                        geom_line(aes(col = Condition))+
                                        facet_wrap(~Channel)+
                                        labs(labs)+
                                        scale_color_manual(values = values ,
                                                           name = colnames(dta)[test])+
                                        scale_fill_manual(values = values ,
                                                          name = colnames(dta)[test])+
                                        ylim(ylim)+
                                        theme_default()
                                
                        }
                        if (testmode == "numer") {
                                numerlist <- list()
                                datalist <- list()
                                for (i in 1:length(levels(dta[,channel]))){
                                        numerlist[[i]] <- filter(dta,dta[,channel] == levels(dta[,channel])[i])
                                        cormat <- matrix(NA,length(uV),boot.num)
                                        for (j in 1:boot.num){
                                                idx <- sample.int(nrow(numerlist[[i]]), nrow(numerlist[[i]]), replace = TRUE) 
                                                cormat[,j] <- melt(cor(numerlist[[i]][idx,test],numerlist[[i]][idx,uV]))[,3]
                                        }
                                        a <- apply(cormat, 1, quantile, probs = boot.intval)
                                        datalist[[i]] <- data.frame(Channel = levels(dta[,channel])[i],
                                                                    frames = frames,
                                                                    Correlation = melt(cor(numerlist[[i]][,test], numerlist[[i]][,uV]))[,3],
                                                                    corQ1 = a[1,], corQ2 = a[2,])
                                }
                                data_plot <- do.call("rbind",datalist)
                                if (order == TRUE){
                                        dta_clustorder <- dcast(data_plot[,1:3],Channel ~ frames,value.var = "Correlation")
                                        hc.pred <- hclust(diss(as.matrix(dta_clustorder[,-1]), "ACF"))
                                        dta_clustorder <- dta_clustorder[hc.pred$order,]
                                        data_plot$Channel <- factor(data_plot$Channel,levels = c(as.character(dta_clustorder$Channel)))
                                }
                                plot <- ggplot(data_plot,
                                               aes(x =frames, y = Correlation))+
                                        geom_ribbon(aes(ymax = corQ2,ymin=corQ1),alpha = ci.alpha,fill= curve.col[1])+
                                        geom_line(col = curve.col[1]) +
                                        labs(labs) +
                                        ylim(ylim) +
                                        theme_bw()+
                                        facet_wrap(~Channel)+
                                        theme_default()
                                
                                #stop("Numeric variable could not adopt bootci ! It will be meaningless.")
                        }
                }
                print(plot)
        }
        if (scalp == T){ # do.call  grid.arrange problem
                message("Label of X-axis will be suppressed.")
                get_legend<-function(myggplot){
                        tmp <- ggplot_gtable(ggplot_build(myggplot))
                        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                        legend <- tmp$grobs[[leg]]
                        return(legend)
                }
                if (is.null(coord.mat)){
                        message("Use default coord.mat")
                        erplay <- rbind(c(NA,NA,NA,NA,"Fp1","Fpz","Fp2",NA,NA,NA,NA),
                                        c(NA,NA,NA,"AF7","AF3","AFz","AF4","AF8",NA,NA,NA),
                                        c(NA,"F7","F5","F3","F1","Fz","F2","F4","F6","F8",NA),
                                        c("FT9","FT7","FC5","FC3","FC1","FCz","FC2","FC4","FC6","FT8","FT10"),
                                        c("T9","T7","C5","C3","C1","CZ","C2","C4","C6","T8","T10"),
                                        c("TP9","TP7","CP5","CP3","CP1","CPz","CP2","CP4","CP6","TP8","TP10"),
                                        c("P9","P7","P5","P3","P1","PZ","P2","P4","P6","P8","P10"),
                                        c(NA,NA,"PO9","PO7","PO3","POz","PO4","PO8","PO10",NA,NA),
                                        c(NA,NA,NA,NA,"O1","Oz","O2",NA,NA,NA,NA),
                                        c(NA,NA,NA,NA,"I1","Iz","I2",NA,NA,NA,NA))
                } else { erplay <- coord.mat } 
                plotlist <- list()
                if (mode == "raw"){
                        data_plot <- dta
                        colnames(data_plot)[channel] <- "Channel"
                        plotlist <- list()
                        if (testmode == "single"){
                                for (i in 1:length(levels(data_plot[,channel]))){
                                        dta_sub <- subset(data_plot,data_plot[,channel] == levels(data_plot[,channel])[i])
                                        plotlist[[i]] <- edaplot(data = dta_sub ,frames = frames,uV = uV,subject=subject,
                                                                 channel=channel)+
                                                geom_line(col = curve.col[1])+
                                                labs(labs(list(title=levels(data_plot[,channel])[i])))+
                                                ylim(ylim)+
                                                theme_default()
                                        names(plotlist)[i] <-levels(data_plot[,channel])[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                                string <- as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        print(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        }
                        if (testmode == "typical"){
                                colnames(data_plot)[test] <- "Condition"
                                values <- curve.col[1]
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                dta_sub <- subset(data_plot,data_plot[,channel] == levels(data_plot[,channel])[1])
                                flegend <- edaplot(data = dta_sub ,frames = frames,uV = uV,subject=subject,
                                                   channel=channel)+
                                        geom_line(aes(col = Condition))+
                                        labs(labs(list(title=levels(data_plot[,channel])[1])))+
                                        ylim(ylim)+
                                        scale_color_manual(values = values, name = colnames(dta)[test])+
                                        theme_default()
                                legend <- get_legend(flegend)
                                for (i in 1:length(levels(data_plot[,channel]))){
                                        dta_sub <- subset(data_plot,data_plot[,channel] == levels(data_plot[,channel])[i])
                                        plotlist[[i]] <- edaplot(data = dta_sub ,frames = frames,uV = uV,subject=subject,
                                                                 channel=channel)+
                                                geom_line(aes(col = Condition))+
                                                labs(labs(list(title=levels(data_plot[,channel])[i])))+
                                                ylim(ylim)+
                                                scale_color_manual(values = values)+
                                                theme_default()+
                                                theme(legend.position = "none")
                                        names(plotlist)[i] <-levels(data_plot[,channel])[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                                string <- as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                printlist[[length(printlist)+1]] <- legend
                                erplaynum[1,ncol(erplaynum)] <- length(printlist) 
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        }
                        if (testmode == "numer"){
                                colnames(data_plot)[test] <- "numer"
                                dta_sub <- subset(data_plot,data_plot[,channel] == levels(data_plot[,channel])[1])
                                flegend <- edaplot(data = dta_sub ,frames = frames,uV = uV,subject=subject,
                                                   channel=channel)+
                                        geom_line(aes(col = numer))+
                                        labs(labs(list(title=levels(data_plot[,channel])[1])))+
                                        ylim(ylim)+
                                        scale_color_gradient(low =curve.col[1], high = curve.col[2],name =  colnames(data)[test])+
                                        theme_default()
                                legend <- get_legend(flegend)
                                for (i in 1:length(levels(data_plot[,channel]))){
                                        dta_sub <- subset(data_plot,data_plot[,channel] == levels(data_plot[,channel])[i])
                                        plotlist[[i]] <- edaplot(data = dta_sub ,frames = frames,uV = uV,subject=subject,
                                                                 channel=channel)+
                                                geom_line(aes(col = numer))+
                                                labs(labs(list(title=levels(data_plot[,channel])[i])))+
                                                ylim(ylim)+
                                                scale_color_gradient(low =curve.col[1], high = curve.col[2],name =  colnames(data)[test])+
                                                theme_default()+
                                                theme(legend.position = "none")
                                        names(plotlist)[i] <-levels(data_plot[,channel])[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                                string <- as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                printlist[[length(printlist)+1]] <- legend
                                erplaynum[1,ncol(erplaynum)] <- length(printlist) 
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        }
                }
                if (mode == "mean") {
                        plotlist <- list()
                        if (testmode == "single"){
                                data_plot <- data_summarize(dta,uV ,summary.var = channel,fun = curve.fun)
                                data_plot <-  suppressMessages(melt(data_plot))
                                colnames(data_plot)[1] <- "Channel"
                                data_plot <- data_plot[order(data_plot$Channel,data_plot$variable),]
                                data_plot$frames <- rep(frames,dim(data_plot)[1]/length(frames))
                                for (i in 1:length(levels(data_plot[,1]))){
                                        dta_sub <- subset(data_plot,data_plot[,1] == levels(data_plot[,1])[i])
                                        plotlist[[i]] <- ggplot(data = dta_sub,aes(x = frames, y = value))+
                                                geom_line(col = curve.col[1])+
                                                theme_bw()+
                                                ylim(ylim)+
                                                labs(labs(list(title=levels(data_plot[,channel])[i])))+
                                                theme_default()
                                        names(plotlist)[i] <-levels(data_plot[,1])[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                                string <- as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        } 
                        if (testmode == "typical") {
                                values <- curve.col[1]
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                data_plot <- data_summarize(dta,uV ,summary.var = c(channel,test),fun = curve.fun)
                                data_plot <-  suppressMessages(melt(data_plot))
                                colnames(data_plot)[c(1,2)] <- c("Channel","Condition")
                                data_plot <- data_plot[order(data_plot$Channel,data_plot$Condition,data_plot$variable),]
                                data_plot$frames <- rep(frames,dim(data_plot)[1]/length(frames))
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                dta_sub <- subset(data_plot,data_plot[,1] == levels(data_plot[,1])[1])
                                flegend <- ggplot(data = dta_sub,aes(x = frames, y = value))+
                                        geom_line(aes(col = Condition))+
                                        labs(labs(list(title=levels(data_plot[,1])[1])))+
                                        ylim(ylim)+
                                        scale_color_manual(values = values, name = colnames(dta)[test])+
                                        theme_default()
                                legend <- get_legend(flegend)
                                for (i in 1:length(levels(data_plot[,1]))){
                                        dta_sub <- subset(data_plot,data_plot[,1] == levels(data_plot[,1])[i])
                                        plotlist[[i]] <- ggplot(data = dta_sub,aes(x = frames, y = value))+
                                                geom_line(aes(col = Condition))+
                                                labs(labs(list(title=levels(data_plot[,1])[i])))+
                                                ylim(ylim)+
                                                scale_color_manual(values = values, name = colnames(dta)[test])+
                                                theme_default()+
                                                theme(legend.position = "none")
                                        names(plotlist)[i] <-levels(data_plot[,1])[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_plot$Channel))] = NA
                                string <- as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_plot$Channel[!toupper(data_plot$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_plot$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                printlist[[length(printlist)+1]] <- legend
                                erplaynum[1,ncol(erplaynum)] <- length(printlist)
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        }
                        if (testmode == "numer") {stop("Numeric variable could not adopt mean curve ! It will be meaningless.")}
                } 
                if (mode == "bootci") {
                        plotlist <- list()
                        if (testmode == "single"){
                                data_ci <- data_ci[order(data_ci$Channel),]
                                data_ci$frames <- rep(frames,dim(data_ci)[1]/length(frames))
                                for (i in 1:length(levels(data_ci$Channel))){
                                        dta_sub <- subset(data_ci,data_ci$Channel == levels(data_ci$Channel)[i])
                                        plotlist[[i]] <- ggplot(data = dta_sub,aes(x = frames, y = FUN))+
                                                geom_ribbon(aes(ymax = Q1,ymin = Q2),alpha = ci.alpha, fill = curve.col[1])+
                                                geom_line(col = curve.col[1])+
                                                ylim(ylim)+
                                                labs(labs(list(title=levels(data_ci[,1])[i])))+
                                                theme_default()
                                        names(plotlist)[i] <-levels(data_ci$Channel)[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_ci$Channel))] = NA
                                string <- as.character(data_ci$Channel[!toupper(data_ci$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_ci$Channel[!toupper(data_ci$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_ci$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        } 
                        if (testmode == "typical") {
                                values <- curve.col[1]
                                for (i in 2:length(levels(as.factor(dta[,test])))){
                                        values <- c(values,curve.col[i])
                                }
                                data_ci <- data_ci[order(data_ci$Channel,data_ci$Condition,data_ci$variable),]
                                data_ci$frames <- rep(frames,dim(data_ci)[1]/length(frames))
                                dta_sub <- subset(data_ci,data_ci[,1] == levels(data_ci[,1])[1])
                                flegend <- ggplot(data = dta_sub,aes(x = frames, y = FUN))+
                                        geom_ribbon(aes(ymax = Q2,ymin=Q1, fill = Condition),alpha = ci.alpha)+
                                        geom_line(aes(col = Condition))+
                                        labs(labs(list(title=levels(data_ci[,1])[1])))+
                                        ylim(ylim)+
                                        scale_color_manual(values = values, name = colnames(dta)[test])+
                                        scale_fill_manual(values = values ,name = colnames(dta)[test])+
                                        theme_default()
                                legend <- get_legend(flegend)
                                for (i in 1:length(levels(data_ci$Channel))){
                                        dta_sub <- subset(data_ci,data_ci[,1] == levels(data_ci[,1])[i])
                                        plotlist[[i]] <- ggplot(data = dta_sub,aes(x = frames, y = FUN))+
                                                geom_ribbon(aes(ymax = Q2,ymin=Q1, fill = Condition),alpha = ci.alpha)+
                                                geom_line(aes(col = Condition))+
                                                labs(labs(list(title=levels(data_ci[,1])[i])))+
                                                ylim(ylim)+
                                                scale_color_manual(values = values, name = colnames(dta)[test])+
                                                scale_fill_manual(values = values ,name = colnames(dta)[test])+
                                                theme_default()+
                                                theme(legend.position = "none")
                                        names(plotlist)[i] <-levels(data_ci$Channel)[i]
                                }
                                erplay[(!toupper(erplay) %in% toupper(data_ci$Channel))] = NA
                                string <- as.character(data_ci$Channel[!toupper(data_ci$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_ci$Channel[!toupper(data_ci$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_ci$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                printlist[[length(printlist)+1]] <- legend
                                erplaynum[1,ncol(erplaynum)] <- length(printlist)
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                        } 
                        if (testmode == "numer"){
                                numerlist <- list()
                                datalist <- list()
                                plotlist <- list()
                                for (i in 1:length(levels(dta[,channel]))){
                                        numerlist[[i]] <- filter(dta,dta[,channel] == levels(dta[,channel])[i])
                                        cormat <- matrix(NA,length(uV),boot.num)
                                        for (j in 1:boot.num){
                                                idx <- sample.int(nrow(numerlist[[i]]), nrow(numerlist[[i]]), replace = TRUE) 
                                                cormat[,j] <- melt(cor(numerlist[[i]][idx,test],numerlist[[i]][idx,uV]))[,3]
                                        }
                                        a <- apply(cormat, 1, quantile, probs = boot.intval)
                                        datalist[[i]] <- data.frame(Channel = levels(dta[,channel])[i],
                                                                    frames = frames,
                                                                    Correlation = melt(cor(numerlist[[i]][,test], numerlist[[i]][,uV]))[,3],
                                                                    corQ1 = a[1,], corQ2 = a[2,])
                                        plotlist[[i]] <- ggplot(datalist[[i]],aes(x =frames, y = Correlation))+
                                                geom_ribbon(aes(ymax = corQ2,ymin=corQ1),alpha = ci.alpha,fill= curve.col[1])+
                                                geom_line(col = curve.col[1]) +
                                                labs(labs) +
                                                ylim(ylim) +
                                                theme_bw()+
                                                facet_wrap(~Channel)+
                                                theme_default()
                                        names(plotlist)[i] <-levels(dta[,channel])[i]
                                }
                                data_ci <- do.call("rbind",datalist)
                                erplay[(!toupper(erplay) %in% toupper(data_ci$Channel))] = NA
                                string <- as.character(data_ci$Channel[!toupper(data_ci$Channel) %in% toupper(erplay)])[!duplicated(as.character(data_ci$Channel[!toupper(data_ci$Channel) %in% toupper(erplay)]))]
                                if (!identical(string, character(0))){
                                        message("The following electrodes will not be plotted :")
                                        message(string)
                                }
                                printlist <- plotlist[which(toupper(names(plotlist)) %in% toupper(erplay))]
                                erplaynum <- matrix(NA,nrow(erplay),ncol(erplay))
                                for (i in 1:nrow(erplay)){
                                        for (j in 1:ncol(erplay)){
                                                if (toupper(erplay[i,j]) %in% toupper(data_ci$Channel)){
                                                        erplaynum[i,j] = which(toupper(names(printlist)) %in% toupper(erplay[i,j]))
                                                }
                                        }
                                }
                                plot <- do.call("grid.arrange", list(grobs=printlist,layout_matrix=erplaynum))
                                #stop("Numeric variable could not adopt bootci ! It will be meaningless.")
                        }
                }
        }
        return(plot)
}
