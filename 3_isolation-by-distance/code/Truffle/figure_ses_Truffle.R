rm(list=ls())
#' --------------------------------------------------------------------------   @Header
#'
#' @title Figure 2 - Panel of Slope with information about SES values for each species
#'
#' @description
#' This R script...
#'
#' @author   Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/02/03
#'
#' --------------------------------------------------------------------------   @VVVVVV
setwd("/Users/brigittetombers/Desktop/DistancesR/R/Truffle")


#'  -------------------------------------------------------------------------   @library
library(ggplot2)
library(gridExtra)
library(cowplot)
library(png)
#'  -------------------------------------------------------------------------   @Parameters
species_names<-c("H_puella")


#'  -------------------------------------------------------------------------   @InitExportDevice

data_names<-list.files(pattern="*Rdata",recursive = FALSE)
lapply(data_names,load,.GlobalEnv)
slope<-list(Slope_SES)

data_pic<-list.files(pattern = "*.png",recursive = FALSE)
pic<-lapply(data_pic,readPNG)
names(pic) <- species_names
#'  -------------------------------------------------------------------------   @formatdata
prep.data <-  function(data){ 
  new_data <- data.frame(windows = data[[1]][1],
             slope_m = apply(do.call(cbind,lapply(data,function(x) x[,2])),1,mean),
             slope_sd = apply(do.call(cbind,lapply(data,function(x) x[,2])),1,sd),
             ses_m = apply(do.call(cbind,lapply(data,function(x) x[,4])),1,mean),
             ses_sd = apply(do.call(cbind,lapply(data,function(x) x[,4])),1,sd),
             signif= rep(NA,length(data[[1]][1])))
  
 # for(i in 1:nrow(new_data)){
 #   if (abs(new_data$ses_m[i])<1.96) new_data$signif[i] <- "nf"
 #   if (abs(new_data$ses_m[i])>1.96) new_data$signif[i] <- "signif"
 # }
  
  return(new_data)
}

data_graph<-lapply(slope,prep.data) 
names(data_graph) <- species_names

#'  -------------------------------------------------------------------------   @DrawFigure
#'  SAVE in 10 per 12
#--- Diplodus graph
        #--- Plot of slope
        plot_H_puella_slope <- ggplot(data_graph[["H_puella"]], aes(x=windows, y=slope_m)) +
          scale_y_continuous(labels=seq(-0.0005,0.0001,0.0001),limits=c(-0.0005,0.0001),
                             breaks = seq(-0.0005,0.0001,0.0001))+
          scale_shape_manual(values=c(16,21))+
          scale_size_manual(values=c(1.2,3))+
          scale_color_manual(values = c("nf" = "#EE00EE32","signif" = "magenta2"))+
          scale_fill_manual( values = c("nf" = "#EE00EE32","signif" = "grey50"))+
          ylab("Slope")+
          xlab(" ")+
          xlim(0,300)+
          geom_errorbar(aes(ymax= slope_m + slope_sd, ymin = slope_m - slope_sd,color=signif,size=signif),size=0.5)+
          geom_point(aes(color=signif,shape=signif,fill=signif,size=signif))+
          theme_bw()+theme(legend.position = "none",
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
         (plot_H_Puella_slope <- ggdraw(plot_H_puella_slope)+draw_image(pic$H_puella,x=0.38,y=0.001,vjust=0.3,scale=0.2))
          
          
        
        #--- Plot of SES    
        plot_H_puella_ses <- ggplot(data_graph[["H_puella"]], aes(x=windows, y=ses_m)) +
          scale_shape_manual(values=c(16,21))+
          scale_size_manual(values=c(1.2,3))+
          scale_color_manual(values = c("nf" = "#EE00EE32","signif" = "magenta2"))+
          scale_fill_manual( values = c("nf" = "#EE00EE32","signif" = "grey50"))+
          ylab("SES")+
          xlab(" ")+
          xlim(0,250)+
          ylim(-6,3)+
          geom_errorbar(aes(ymax= ses_m + ses_sd, ymin = ses_m - ses_sd,color=signif,size=signif),size=0.5)+
          geom_hline(yintercept=0,  color = "black")+ 
          geom_hline(yintercept=1.96,linetype="dashed", color = "black")+
          geom_hline(yintercept=-1.96, linetype="dashed", color = "black")+
          geom_point(aes(color=signif,shape=signif,fill=signif,size=signif))+
          theme_bw()+theme(legend.position = "none",
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        (plot_H_puella_ses <- ggdraw(plot_H_puella_ses)+draw_image(pic$H_puella,x=0.38,y=0.001,vjust=0.3,scale=0.2))
        
       
          
        ggsave(filename = "plot_H_puella_ses_truffle_2.5k_1k.pdf", plot = plot_H_puella_ses)
        ggsave(filename = "plot_H_puella_slope_truffle_2.5k_1k.pdf", plot = plot_H_puella_slope) 
