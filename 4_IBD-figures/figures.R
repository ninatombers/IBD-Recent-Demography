rm(list=ls())

library("tidyverse")
library(ggplot2)
library(plyr)
library(dplyr)


setwd("/Users/brigittetombers/Desktop/Project/Forschungsprojekt_Figures")
 
#----------------NUMBER OF SEGMENTS---------------------------------------------------------------
segments_truffle <- read_tsv("truffle_1k_0.5k_nofilt_nocontig.segments.tsv") 
segments_tpbwt <- read_csv("300_0.01.csv") %>%
  mutate(length_bp=(end_bp-start_bp)) 

number_segments <- data.frame(LG = c("LG01","LG02","LG03","LG04","LG05","LG06","LG07","LG08","LG09","LG10","LG11","LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21","LG22","LG23","LG24"))

#second cloumn Truffle number of segment
for (i in 1:24) {
  
  number_segments[i,2] <- sum(segments_truffle$CHROM == number_segments[i,1])
}

#second cloumn TPBWT number of segment
for (i in 1:24) {
  
  number_segments[i,3] <- sum(segments_tpbwt$chromosome == number_segments[i,1])
}

#name coloumns 
colnames(number_segments) <- c("LG","nsegmt_truffle","nsegmt_TPBWT")

##number of segments over all
sum(number_segments$nsegmt_truffle)
sum(number_segments$nsegmt_TPBWT)

# number of IBD1 /IBD2
sum(segments_truffle$TYPE == "IBD1")
sum(segments_truffle$TYPE == "IBD2")

# average segment length (TRUFFLE Mbp)
sum(segments_tpbwt$length_bp)/sum(number_segments$nsegmt_TPBWT)
sum(segments_truffle$LENGTH)/sum(number_segments$nsegmt_truffle)*10^6

#Graph mit Segmenten pro LG
nSeg1 <- data.frame(LG = number_segments[, 1])
nSeg1[, 2] <- data.frame(nSegments = number_segments[, 2])
nSeg2 <- data.frame(LG = number_segments[, 1])
nSeg2[, 2] <- data.frame(nSegments = number_segments[, 3])


nSeg1$method <- 'TRUFFLE'
nSeg2$method <- 'TPBWT'

nSeg <- rbind(nSeg1, nSeg2)

p <- ggplot(nSeg, mapping = aes(x = LG, y = nSegments, fill = method)) +
  geom_bar(stat="identity", color="black", position=position_dodge(width=0.5))+
  scale_fill_manual(values=c(rgb(0.2,0.8,1, 0.4), rgb(1,0.6,0, 0.4)))  + 
  labs(x = "LinkageGroup") +
  theme_bw(base_size = 12, base_family = "") + 
  labs(title = "Number of Segments per LG")+
  theme(plot.title = element_text(face="bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
 
#-------------------------LENGTH OF SEGMENTS------------------------

etapa1 <- data.frame(Segmlength = log((segments_truffle$LENGTH)*10^6))
etapa2 <- data.frame(Segmlength = log(segments_tpbwt$length_bp))

etapa1$method <- 'TRUFFLE'
etapa2$method <- 'TPBWT'

combo <- rbind(etapa1, etapa2)

plot <- ggplot(combo, mapping = aes(x = Segmlength, fill = method)) + 
      geom_histogram(color="black", position = "identity",alpha = 0.4, bins = 100) +
      scale_fill_manual(values=c(rgb(0.2,0.8,1,0.4), rgb(1,0.6,0,0.4))) + 
      theme(legend.position="bottom") + 
      labs(x = "log(segment_lentgh)[Mbp]") +
      theme_bw(base_size = 12, base_family = "") + 
      labs(title = "Histogram Segment Length")+
      theme(plot.title = element_text(face="bold")) +
      theme(plot.title = element_text(hjust = 0.5))


