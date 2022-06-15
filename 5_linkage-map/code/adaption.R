##This file put the output of TPBWT in the format as truffle so the conversion 
##script to cM can be used
## by: Nina Tombers: 01/03/2022
rm(list=ls())

library("tidyverse")
library("superheat")

setwd("/Users/brigittetombers/Desktop/Project/5.Linkage-map/TPBWT")

#read in data and caculate length_bp
#data_fastIBD <- read_delim("300_0.2.csv", delim = ",") %>%
#data_fastIBD <- read_delim("300_0.01.csv", delim = ",") %>%
#data_fastIBD <- read_delim("TPBWT_0.1.csv", delim = ",") %>%
data_fastIBD <- read_delim("300_0.01_noContig.csv", delim = ",") %>%
  #data_fastIBD <- read_delim("300_0.15.csv", delim = ",") %>%
  #data_fastIBD <- read_delim("250_0.2.csv", delim = ",") %>%
  #data_fastIBD <- read_delim("100_0.3.csv", delim = ",") %>%
  mutate(length_bp=(end_bp-start_bp)) 

# make .IBD file------------------------------------------

#connect id1 and id2 from phasedIBD with list "samplenames" sample names

samplenames <- read_delim("samplenames.csv", delim = ";")
samplenamesid2 <- rename(samplenames, id2 = id1, ID2 = ID1)

data_fastIBD <- left_join(data_fastIBD, samplenames) 
data_fastIBD <- left_join(data_fastIBD, samplenamesid2)

#calculate IBD-values from bp-length, by summing length over all with same ID1&ID2, then using this (IBD2 + 0.5·IBD1)/(IBD0 +IBD1 +IBD2) in which each value in data is IBD1

ibd_pI <- select(data_fastIBD, ID1, ID2, length_bp) %>%
  arrange(ID1, ID2)

ibd_new <- ibd_pI %>% 
  group_by(ID1, ID2) %>% 
  summarise_all(funs(sum)) %>%
  mutate(IBD1=((length_bp*0.5)/559649677))



ibd_new$IBD1 <- as.numeric(ibd_new$IBD1)
ibd_new <- select(ibd_new, ID1, ID2, IBD1)

TPBWT_ibd <- ibd_new %>% 
  mutate(IBD0=(1-IBD1), IBD2="0", SEX="00")

# make .segments file------------------------------------------

TPBWT_seg <- data_fastIBD %>% 
  mutate(TYPE="IBD1", ID1=id1, ID2=id2, CHROM=chromosome, VARSTART=start_bp, VAREND=end_bp, POS=start_bp, LENGTH=(end_bp-start_bp), NMARKERS=end-start) %>%
  select(TYPE,ID1,ID2,CHROM,VARSTART,VAREND,POS,LENGTH,NMARKERS)


write.table(TPBWT_ibd, file = "300_0.01.ibd", sep = ",", quote = FALSE, row.names = F)
write.table(TPBWT_seg, file = "300_0.01.segments", sep = ",", quote = FALSE, row.names = F)
