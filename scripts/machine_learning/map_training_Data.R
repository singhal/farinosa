rm(list = ls())

library(ggplot2)
library(sf)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/")

d1 = list.files("machine_learning/inat_observations/", 
                full.names = T)
d2 = lapply(d1, "read.csv")
cols = table(unlist(lapply(d2, names)))
cols = names(cols[ cols == length(d1)])
d3 = lapply(d2, function(x) x[, cols])
d3 = do.call("rbind", d3)

# load in training data
x = read.csv("machine_learning/Encelia.training_validation.csv")
# lead in accuracy data
y = read.csv("machine_learning/accuracy_set.csv")

x1 = x %>% select(file = file, flowers = s) %>% 
  rbind(y %>% select(file, flowers))
x2 = x1[x1$flowers %in% c("B", "Y", "N"), ]
x2$id = gsub(".jpg", "", x2$file)

x3 = cbind(x2, d3[match(x2$id, d3$id), c("latitude", "longitude", "observed_on")])

world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
   geom_sf() +
   geom_point(data = x3 %>% filter(flowers != "N"), 
              aes(x = longitude, 
                  y = latitude,
                  color = flowers), 
              shape = 16, alpha = 1) +
  coord_sf(xlim = c(-120, -108), ylim = c(23, 37)) +
  scale_color_manual(values = c("black", "yellow")) +
  theme_minimal()
   
ggplot(data = world) +
    geom_sf() +
    geom_point(data = x3, 
               aes(x = longitude, 
                   y = latitude), 
               shape = 16, alpha = 0.2) +
    coord_sf(xlim = c(-120, -108), ylim = c(23, 37)) +
  theme_minimal()

m = x3 %>% 
  mutate(month = gsub("-\\S\\S", "", 
                      gsub("^\\S\\S\\S\\S-", "", x3$observed_on))) %>%
  select(month, flowers) %>% group_by(month) %>% 
  summarize(Y_flowering = sum(flowers == "Y") / n(), 
            B_flowering = sum(flowers == "B") / n())
ggplot(m) +
  geom_point(aes(x = month, y = Y_flowering), pch = 21, fill = "yellow") +
  geom_point(aes(x = month, y = B_flowering), pch = 21, fill = "brown") +
  theme_bw()
