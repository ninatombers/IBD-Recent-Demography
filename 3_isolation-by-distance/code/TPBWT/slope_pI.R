require(ggplot2)
require(grid)
require(gridExtra)
require(reshape2)

setwd("/Users/brigittetombers/Desktop/DistancesR/R/PhasedIBD")

colnames(final_pI)[1] <- "INDV1"
colnames(final_pI)[2] <- "INDV2"
colnames(final_pI)[3] <- "DIST"
colnames(final_pI)[4] <- "Ajk"

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)

source("/Users/brigittetombers/Desktop/DistancesR/R/PhasedIBD/Slope.boot.SES_pI.R")

#min <- nrow(subset(final, final$DIST <5))

###Slope_SES <- Slope.boot.SES(final, window = c(seq(from = 5,to = max(final$DIST),by=5)), min = nrow(subset(final, final$DIST <5)), repnull = 1000, repboot = 1000,core = 2)

Slope_SES <- Slope.boot.SES(final_pI, window = c(seq(from = 5,to = max(final_pI$DIST),by=5)), min = nrow(subset(final_pI, final_pI$DIST <5)), repnull = 100, repboot = 100,core = 2)
#Slope_SES <- Slope.boot.SES(final, window = c(seq(5,920,by=5)),min = min,repnull = 100,repboot = 100,core = 4)

Slope_H_Puella <- data.frame(aaply(laply(Slope_SES, as.matrix), c(2, 3), mean),aaply(laply(Slope_SES, as.matrix), c(2, 3), sd)[,c(2,4)])
colnames(Slope_H_Puella)[c(5,6)]<- c("slope_sd","SES_sd")
                       

### Save the results in a file
write.table(Slope_H_Puella, "Slope_H_Puella.txt", quote=FALSE, row.names=FALSE)

### Save the results in a Rdatasave.image("Slope_H_Puella.RData")
save.image("Slope_H_puella_pI.RData")

