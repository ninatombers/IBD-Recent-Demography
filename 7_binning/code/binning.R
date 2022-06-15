rm(list=ls())

library(tidyverse)
library(ggplot2)
library(OneR)


setwd("/Users/brigittetombers/Desktop/Project/6.Binning")

# excluding regions >95% IBD score
#truffle_seg_filt <- read_tsv("truffle_1k0.5k_95_segments_filt.tsv") %>%
# all regions
  truffle_seg_filt <- read_tsv("truffle_1k_0.5k.segments_filt.tsv") %>%
  #seperate PAIR into ID1 & ID2
mutate(ID1=substr(PAIR,1,11)) %>%
mutate (ID2=substr(PAIR,13,23))

#truffle_seg_filt <-filter(truffle_seg_filt,truffle_seg_filt$CHROM!="LG08")

pairs <- unique(unique(truffle_seg_filt[,c("ID1", "ID2")]))
#----------------bin distance-----------------------

#load geographic distance
geo <- read_delim("geographical-distance.csv", delim = ";")
lgeo <- pivot_longer(geo, !ID, names_to = "pair", values_to = "distance") %>%
  select(pair, ID, distance) %>%
  arrange(pair, ID) %>%
  rename(ID1 = pair, ID2 = ID) %>%
  drop_na()
#add geographic distance to table
truffle_con_geo <- left_join(lgeo, truffle_seg_filt)

#count number of pairs per distance
library(plyr)
freq_pairs <- count(lgeo$distance)
colnames(freq_pairs) <- c("distance", "freq_pair")

#------------------bin length-----------------------------

#count = number of segments per distance

# bin into different cM lengths map1
#filter 
#bin1
bin1 <-   truffle_con_geo%>%
  filter(length_cM_m1 >= 0.2 & length_cM_m1 < 0.8  )
bin1 <- full_join(bin1, pairs)

sum1 <- count(bin1,vars = c("PAIR","distance")) %>% 
  group_by(distance) %>%
  summarise_each(funs(mean(., na.rm=T), sd, n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))), freq)

#count number of segments in this bin for each distance
count1 <- count(bin1$distance)
colnames(count1) <- c("distance", "freq_seg")
count1 <- left_join(count1, sum1)


#bin2
bin2 <-   truffle_con_geo%>%
  filter(length_cM_m1 >= 0.8 & length_cM_m1 < 1.5 )
bin2 <- full_join(bin2, pairs)

sum2 <- count(bin2,vars = c("PAIR","distance")) %>% 
  group_by(distance) %>%
  summarise_each(funs(mean(., na.rm=T), sd, n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))), freq)

count2 <- count(bin2$distance)
colnames(count2) <- c("distance", "freq_seg")
count2 <- left_join(count2, sum2)

#bin3
bin3 <-   truffle_con_geo%>%
  filter(length_cM_m1 >= 1.5 & length_cM_m1 < 2.9 )
bin3 <- full_join(bin3, pairs)

sum3 <- count(bin3,vars = c("PAIR","distance")) %>% 
  group_by(distance) %>%
  summarise_each(funs(mean(., na.rm=T), sd, n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))), freq)

count3 <- count(bin3$distance)
colnames(count3) <- c("distance", "freq_seg")
count3 <- left_join(count3, sum3)

#bin4
bin4 <-   truffle_con_geo%>%
  filter(length_cM_m1 >= 2.9)
bin4 <- full_join(bin4, pairs)

sum4 <- count(bin4,vars = c("PAIR","distance")) %>% 
  group_by(distance) %>%
  summarise_each(funs(mean(., na.rm=T), sd, n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))), freq)

count4 <- count(bin4$distance)
count4 <- count(bin4, vars = "distance")
colnames(count4) <- c("distance", "freq_seg")
count4 <- left_join(count4, sum4)


# join number of segments and pairs for each cM bin and add collum with bin
b1 <-left_join(count1, freq_pairs) %>%
    mutate(bin="0.2 - 0.8 cM")
b2 <-left_join(count2, freq_pairs) %>%
  mutate(bin="0.8 - 1.5 cM")
b3 <-left_join(count3, freq_pairs) %>%
  mutate(bin="1.5.0 - 2.9 cM")
b4 <-left_join(count4, freq_pairs) %>%
  mutate(bin="2.9+cM")


#merge into one big file with all bins, number segments & number pairs
library(dplyr)
bins <- rbind(b1,b2,b3,b4)

dist <- data.frame(unique(lgeo$distance))
colnames(dist) <- c("distance")

#control:
#new <- mutate(bins, y=(freq_seg/freq_pair)) 
sum(bins$freq_seg,na.rm = FALSE)

model1 <- lm((freq_seg/freq_pair) ~ distance, data = b1)
model2 <- lm((freq_seg/freq_pair) ~ distance, data = b2)
model3 <- lm((freq_seg/freq_pair) ~ distance, data = b3)
model4 <- lm((freq_seg/freq_pair) ~ distance, data = b4)

#-------------------graph-------------------------

library("gridExtra")
library("cowplot")
require(scales)

(distance_plot <-
   ggplot(bins, aes(x = distance, y = (freq_seg/freq_pair), colour=bin, fill=bin)) +
   geom_point() + 
   geom_smooth(method = 'lm', se=TRUE) +
   theme_bw() + 
    xlab("Distance [km]") +
    ylab("log(nSegments per pair)") +
    ggtitle("Number of segments over Distance in cM_truffle1k_0.5k") +
    geom_errorbar(aes(ymin=(freq_seg/freq_pair)-se, ymax=(freq_seg/freq_pair)+se), width=.2,
                  position=position_dodge(.9)) +
     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="lr")
)


ggsave(filename = "binned_graph_allpairs_zeros_1_0.5>3,5.pdf", plot = distance_plot)

