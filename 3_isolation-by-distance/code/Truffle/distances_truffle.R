rm(list=ls())

library(tidyverse)
library(superheat)


# set directory (with filtering or without)


#with filtering----------------------------------------------------------------#
setwd("/Users/brigittetombers/Desktop/Project/3.DistancesR/R/Truffle")

#data_truffle <- read_delim("truffle_0.5k_0.25k.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_1k_0.5k.csv", delim = ";") %>%
#data_truffle <- read_delim("Truffle_2.5k_1k.csv", delim = ";") %>%
#data_truffle <- read_delim("Truffle_5k_2.5k.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_10k_5k.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_15k_7.5k.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_20k_10k.csv", delim = ";") %>%
#mutate(IBD=((IBD2+0.5*IBD1)/(IBD0+IBD1+IBD2)))

#without filtering -----------------------------------------------------------#

#setwd("/Users/brigittetombers/Desktop/DistancesR/R/Truffle")

# read in truffle data and calculate IBD

data_truffle <- read_tsv("truffle_1k_0.5k_95.ibd.tsv") %>%
#data_truffle <- read_tsv("truffle_1k_0.5k_nofilt_nocontig.ibd.tsv") %>%
  
#data_truffle <- read_delim("Truffle_2.5k_1k_nofilt.csv", delim = ";") %>%
#data_truffle <- read_delim("Truffle_5k_2.5k_nofilt.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_10k_5k_nofilt.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_15k_7.5k_nofilt.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_20k_10k_nofilt.csv", delim = ";") %>%
#data_truffle <- read_delim("truffle_25k_20k_nofilt.csv", delim = ";") %>%
#data_truffle <- read_tsv("truffle_1k_0.5k_nofilt.conv_summary Kopie.tsv") %>%
mutate(IBD=((IBD2+0.5*IBD1)/(IBD0+IBD1+IBD2)))

#-----------------------------------------------------------------------------#

ibd <- select(data_truffle, ID1, ID2, IBD) %>%
  arrange(ID1, ID2)

sampleIDs <- union(data_truffle$ID1, data_truffle$ID2)


# read in geographic distance
geo <- read_delim("geographical-distance.csv", delim = ";")
lgeo <- pivot_longer(geo, !ID, names_to = "pair", values_to = "distance") %>%
  select(pair, ID, distance) %>%
  #filter(distance <= 2.5) %>%
  arrange(pair, ID) %>%
  rename(ID1 = pair, ID2 = ID) %>%
  drop_na()





# create IBD matrix
m <- pivot_wider(ibd, names_from = ID1, values_from = IBD) %>%
  arrange(ID2)

m <- m %>% column_to_rownames(var = 'ID2')


#superheat(m, left.label.text.size = 3, bottom.label.text.size = 1.5, heat.na.col = "white", title = "Superheat",row.title = "ID1",column.title = "ID2")
#add scale = T for  (to mean 0 and standard deviation 1)
#superheat(m, left.label.size = 0.4, bottom.label.size = 0.1, scale = T, heat.na.col = "white", row.dendrogram = TRUE)

# join tables
final <- left_join(lgeo, ibd)

# linear regression model
model <- lm(IBD ~ distance, final)

s <- summary(model)
s$r.squared 

# plot


(distance_plot <-
   ggplot(final, aes(x = distance, y = IBD)) +
   geom_point() +
  geom_smooth(method = 'lm', se=FALSE) +
   theme_minimal() + 
    labs(title = "Isolation by Distance - Truffle 1k/0.5k - 95th Percentile")+
    xlab("Distance (km)") + ylab("IBD score") +
    theme_bw(base_size = 12, base_family = "") +
    theme(plot.title = element_text(hjust = 0.5))
    )

ggsave(filename = "distance_cor_1k_0.5k_nofilter_noContig.pdf", plot = distance_plot)




