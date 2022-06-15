rm(list=ls())

library("tidyverse")
library("superheat")

# set directory (with filtering or without)


setwd("/Users/brigittetombers/Desktop/Project/3.DistancesR/R/Relatedness")

data_truffle <- read.csv("out.relatedness_1k05k.csv")%>%
  arrange(INDV1, INDV2)


#-----------------------------------------------------------------------------#

sampleIDs <- union(data_truffle$INDV1, data_truffle$INDV2)


# read in geographic distance
geo <- read_delim("geographical-distance.csv", delim = ";")
lgeo <- pivot_longer(geo, !ID, names_to = "pair", values_to = "distance") %>%
  select(pair, ID, distance) %>%
  #filter(distance <= 2.5) %>%
  arrange(pair, ID) %>%
  rename(INDV1 = pair, INDV2 = ID) %>%
  drop_na()

# create IBD matrix
m <- pivot_wider(data_truffle, names_from = INDV1, values_from = RELATEDNESS_AJK) %>%
  arrange(INDV2)

m <- m %>% column_to_rownames(var = 'INDV2')


superheat(m, left.label.text.size = 3, bottom.label.text.size = 1.5, heat.na.col = "white", title = "Superheat",row.title = "INDV1",column.title = "INDV2")
#add scale = T for  (to mean 0 and standard deviation 1)
#superheat(m, left.label.size = 0.4, bottom.label.size = 0.1, scale = T, heat.na.col = "white", row.dendrogram = TRUE)

# join tables
final <- left_join(lgeo, data_truffle)

# linear regression model
model <- lm(final$RELATEDNESS_AJK ~ final$distance, data_truffle=final)

s <- summary(model)
s$r.squared 

# plot


(distance_plot <-
   ggplot(final, aes(x = distance, y = RELATEDNESS_AJK)) +
   geom_point() +
  geom_smooth(method = 'lm', se=FALSE) +
   theme_minimal()
    )

ggsave(filename = "distance_cor_1k_0.5k_nofilter_relatedness.pdf", plot = distance_plot)




