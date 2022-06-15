rm(list=ls())

library("tidyverse")
library("superheat")

setwd("/Users/brigittetombers/Desktop/Project/3.DistancesR/R/PhasedIBD")

#read in data and caculate length_bp
#data_fastIBD <- read_delim("300_0.2.csv", delim = ",") %>%
#data_fastIBD <- read_delim("300_0.01.csv", delim = ",") %>%
#data_fastIBD <- read_delim("TPBWT_0.1.csv", delim = ",") %>%
data_fastIBD <- read_delim("300_0.01_noContig.csv", delim = ",") %>%
#data_fastIBD <- read_delim("300_0.15.csv", delim = ",") %>%
#data_fastIBD <- read_delim("250_0.2.csv", delim = ",") %>%
#data_fastIBD <- read_delim("100_0.3.csv", delim = ",") %>%
 mutate(length_bp=(end_bp-start_bp)) 

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
  mutate(IBD=((length_bp*0.5)/559649677))



ibd_new$IBD <- as.numeric(ibd_new$IBD)
ibd_new <- select(ibd_new, ID1, ID2, IBD)



# create IBD matrix
m_pI <- pivot_wider(ibd_new, names_from = ID1, values_from = IBD) %>%
  arrange(ID2)
m_pI <- m_pI %>% column_to_rownames(var = 'ID2')

 
## create heatmap
superheat(m_pI, left.label.text.size = 3, bottom.label.text.size = 1.5, heat.na.col = "white", title = "Superheat",row.title = "ID1",column.title = "ID2")

# read in geographic distance
geo <- read_delim("geographical-distance.csv", delim = ";")
lgeo <- pivot_longer(geo, !ID, names_to = "pair", values_to = "distance") %>%
  select(pair, ID, distance) %>%
  arrange(pair, ID) %>%
  rename(ID1 = pair, ID2 = ID) %>%
  drop_na()


# join tables
final_pI <- left_join(lgeo, ibd_new)

# linear regression model
model <- lm(IBD~distance, final_pI)
#model <- lm(IBD~distance, data_fastIBD=final_pI)

s <- summary(model)
s$r.squared 

# plot


(distance_plot <-
    ggplot(final_pI, aes(x = distance, y = IBD)) +
    geom_point() +
    geom_smooth(method = 'lm', se=FALSE) + 
    theme_minimal()+ 
    labs(title = "Isolation by Distance - TPBWT 300/0.01")+
    xlab("Distance (km)") + ylab("IBD score") +
    theme_bw(base_size = 12, base_family = "") +
    theme(plot.title = element_text(hjust = 0.5))
)




ggsave(filename = "distance_cor_phasedIBD_TPBWT_01.pdf", plot = distance_plot)
