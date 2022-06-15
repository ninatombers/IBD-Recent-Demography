rm(list=ls())

setwd("/Users/brigittetombers/Desktop/Project/9.GenePop")

library("tidyverse")

data_microsat <- read.table("puella_FINAL_GRA.txt", sep=";")

model <- lm(data_microsat$V2 ~ data_microsat$V1,data_microsat)

s <- summary(model)
s$r.squared 

(distance_plot <-
    ggplot(data_microsat, aes(x = V1, y = V2)) +
    geom_point() +
    geom_smooth(method = 'lm', se=FALSE) +
    labs(x = "Distance in km", y = 'Relatedness') +
    theme_bw(base_size = 12, base_family = "") + 
    labs(title = "Isolation by Distance - Mikrosat")+
    theme_minimal()
)
