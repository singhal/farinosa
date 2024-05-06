rm(list = ls())

library(tidyverse)
library(fields)
library(rnaturalearth)
library(scatterpie)
library(patchwork)
library(magick)
library(cowplot)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/")

dd = read.csv("metadata/black_yellow.cleaned11Sept23.csv")
world_map <- ne_countries(scale = "medium", returnclass = "sf")

# https://stackoverflow.com/questions/28672399/spatial-clustering-in-r-simple-example
coors = dd %>% dplyr::select(longitude, latitude)
#distance matrix
dist.in.km.matrix <- rdist.earth(as.data.frame(coors), miles = F, R=6371)
#clustering
fit <- hclust(as.dist(dist.in.km.matrix), method = "single")
# group by 1 vs 20 km
dd = dd %>% mutate(cluster1 = cutree(fit, h = 1), 
                   cluster2 = cutree(fit, h = 20))
length(unique(dd$cluster1))
length(unique(dd$cluster2))

dd20 = dd %>% group_by(cluster2) %>% summarize(yellow=sum(yellow, na.rm=T), 
                                       black = sum(black, na.rm=T), 
                                       latitude = mean(latitude), longitude = mean(longitude)) %>% 
  ungroup() %>% mutate(total = yellow + black) %>%
  mutate(percent_black = black / total)
dd1 = dd %>% group_by(cluster1) %>% summarize(yellow=sum(yellow, na.rm=T), 
                                                         black = sum(black, na.rm=T),
                                              latitude = mean(latitude), longitude = mean(longitude)) %>%
  ungroup() %>% mutate(total = yellow + black) %>%
  mutate(percent_black = black / total)

### zoomed out
inset1x = c(-115.5, -114.8)
inset1y = c(31, 33)

inset2x = c(-117.8, -109.2)
inset2y = c(22.6, 37)
  
A = ggplot() +
  geom_sf(data = world_map, fill = "gray90", color = "gray30") +
  geom_rect(aes(xmin = min(inset1x), xmax = max(inset1x),
                ymin = min(inset1y), ymax = max(inset1y)),
            alpha = 0.5, fill = "gray30") + 
  geom_scatterpie(aes(x=longitude, y=latitude, r=0.15), data=dd20, cols=c("yellow", "black")) +
  scale_fill_manual(values = c(black = "#482624", yellow = "#F5B203")) +
  coord_sf(xlim = inset2x, ylim = inset2y ) +
  theme_void() + theme(legend.position = "none", 
                       panel.border = element_rect(colour = "black", fill=NA, size=0.5))

### zoomed in
B = ggplot() +
  geom_sf(data = world_map, fill = "gray90", color = "gray30") +
  geom_scatterpie(aes(x=longitude, y=latitude, r=0.03), data=dd1, cols=c("yellow", "black")) +
  scale_fill_manual(values = c(black = "#482624", yellow = "#F5B203")) +
  coord_sf(xlim = inset1x, ylim = inset1y) +
  theme_void() + theme(legend.position = "none",
                       panel.border = element_rect(colour = "black", fill=NA, size=0.5))

C = ggdraw() + draw_image("metadata/field.jpg", scale = 0.97)
D = ggdraw() + draw_image("metadata/side_by_side.jpg", scale = 0.97)

inset = ggplot() +
  geom_sf(data = world_map, fill = "gray90", color = "gray90") +
  geom_rect(aes(xmin = min(inset2x), xmax = max(inset2x),
                ymin = min(inset2y), ymax = max(inset2y)),
            col = "black", fill = NA, linewidth = 0.2) +
  coord_sf(xlim = c(-130, -75), ylim = c(10, 80)) +
  theme_void() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                       panel.background = element_rect(fill = "white"))

plt = A + B + ( C / D) + plot_layout(widths = c(2, 1, 2))
save_plot("figures/Figure1_main.pdf", plt, base_height = 4, base_width = 6.5)
save_plot("figures/Figure1_inset.pdf", inset)

