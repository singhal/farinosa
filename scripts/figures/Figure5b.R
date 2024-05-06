rm(list = ls())

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(raster)
library(pcaMethods)
library("rnaturalearth")
library("rnaturalearthdata")

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/")

world = st_read("../../spatial_rasters/ne_50m_land/ne_50m_land.shp")
rivers = st_read("../../spatial_rasters/ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")


# world <- ne_countries(scale = "medium", returnclass = "sf")

d1 = read.csv("machine_learning/accuracy_set.csv") %>%
  filter(!flowers %in% c("U", "M"))

d2 = read.csv("machine_learning/categorization/all_images.lr_0.01.csv")
d2[d2$prediction == 0, "prediction"] = "B"
d2[d2$prediction == 1, "prediction"] = "N"
d2[d2$prediction == 2, "prediction"] = "Y"

# d2 = d2 %>% filter(prediction == "B" & X0 > 0.7 | prediction == "Y" & X2 > 0.7)

d3 = read.csv("machine_learning/Encelia.training_validation.csv")

dd = d1 %>% dplyr::select(file, flowers) %>% 
  rbind(d2 %>% dplyr::select(file, flowers = prediction)) %>%
  rbind(d3 %>% dplyr::select(file, flowers = s))

# lat long data
ll = list.files("machine_learning/inat_observations/", 
                full.names = T)
ll1 = lapply(ll, "read.csv")
cols = table(unlist(lapply(ll1, names)))
cols = names(cols[ cols == length(ll)])
ll2 = lapply(ll1, function(x) x[, cols])
ll2 = do.call("rbind", ll2)
ll3 = ll2 %>% filter(duplicated(id) == FALSE)

dd$id = as.numeric(gsub(".jpg", "", dd$file))
x3 = dd %>% left_join(ll3 %>% dplyr::select(id, latitude, longitude, observed_on))
x4 = x3 %>% filter(latitude > 23)

pts <- st_as_sf(x4 %>% filter(flowers %in% c("B", "Y")),
                coords = c("longitude", "latitude"),
                crs = st_crs(world))

# create a grid covering all points, cell size 0.2 units,
# add unique cell_id
grid_sf <- st_make_grid(pts, cellsize = 0.2) |> 
  st_sf(geometry = _) |>
  mutate(cell_id = row_number(), .before = 1)

# spatial join to metch cell_id-s to points,
grid_sf2 = st_join(pts, grid_sf) %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarize(total = n(), num_black = sum(flowers == "B"),
         per_black = num_black / total) %>%
  ungroup()

# join grid to shares for ploting
grid_sf3 <- right_join(grid_sf, grid_sf2) 


# range map
ll4 = ll3 %>% filter(complete.cases(longitude)) %>% 
  filter(longitude < -108, latitude < 38, latitude > 20) %>%
  dplyr::select(longitude, latitude) %>% summarize_all(round, digits = 3) %>%
  distinct() %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(sf::st_crs(world))
plot(ll4)
all_r1 = alphahull::ashape(st_coordinates(ll4)[, 1:2], alpha = 25)
plot(all_r1)
all_r2 = hull2spatial::ashape2poly(all_r1)
all_r3 = st_as_sf(all_r2) %>% 
  st_set_crs(sf::st_crs(world)) %>% 
  st_intersection(world)

map = ggplot() +
  geom_sf(data = world, fill = "gray95") +
  geom_sf(data = all_r3, fill = "gray80",col = NA) +
  geom_sf(data = rivers, col = "blue4") +
  # draw full grid
  geom_sf(data = grid_sf3 %>% filter(total > 1), aes(fill = per_black), col = NA) +
  coord_sf(xlim = c(-120, -108), ylim = c(23, 37)) +
  theme_minimal() +
  labs(fill = "% brown") +
  scale_fill_gradient(low = "#F5B203", high = "#482624") +
  theme(axis.text = element_blank(), panel.grid = element_blank())
cowplot::save_plot("figures/Figure5.png", map, base_height = 5,
                   base_width = 4)

# ggplot() +
#   geom_sf(data = world, fill = "gray95") +
#   geom_sf(data = all_r3, fill = "gray80",col = NA) +
#   geom_sf(data = rivers, col = "blue4") +
#     # draw full grid
#     geom_sf(data = grid_sf3 %>% filter(num_black > 2, per_black > 0.1), fill = "#482624") +
#     coord_sf(xlim = c(-120, -108), ylim = c(23, 37)) +
#     theme_minimal() +
#     theme(axis.text = element_blank(), panel.grid = element_blank())

# m = x4 %>% 
#   mutate(month = gsub("-\\S\\S", "", 
#                       gsub("^\\S\\S\\S\\S-", "", observed_on))) %>%
#   dplyr::select(month, flowers) %>% group_by(month) %>% 
#   summarize(Y_flowering = sum(flowers == "Y") / n(), 
#             B_flowering = sum(flowers == "B") / n())
# ggplot(m) +
#   geom_point(aes(x = month, y = Y_flowering), pch = 21, fill = "yellow") +
#   geom_point(aes(x = month, y = B_flowering), pch = 21, fill = "brown") +
#   theme_bw()

##########
# climate
#########

climr = list.files("../../spatial_rasters/wc2.0_30s_bio/", pattern=".tif", full.names = T)
climr = stack(lapply(climr, raster))

pts = x4[, c("longitude", "latitude")]
clim = raster::extract(climr, pts)

clim2 = clim[complete.cases(clim), ]
pts2 = x4[complete.cases(clim),]

cc = cbind(clim2, pts2)
xpca = pca(clim2, scale="uv", center = T, nPcs = 10)

summary(xpca)
sort(abs(xpca@loadings[, "PC1"]))
sort(abs(xpca@loadings[, "PC2"]))

pts3 = cbind(xpca@scores[, 1:10], pts2, clim2)

# focus just on the brown flowers at higher latitudes
pts4 = pts3 %>% filter(latitude > 32)

xx = pts4 %>% filter(flowers != "N") %>% group_by(flowers) %>%
  dplyr::select(grep("bio", names(pts4))) %>%
  summarize_all(mean)

ggplot() +
  geom_sf(data = world, fill = "gray95") +
  geom_point(data = pts4 %>% filter(flowers != "N"), aes(longitude, latitude, col = flowers)) +
  coord_sf(xlim = c(-120, -108), ylim = c(32, 37)) +
  scale_color_manual(values = c("#482624", "#F5B203")) +
  theme_classic()

ggplot() +
  geom_point(data = pts4 %>% filter(flowers != "N"), aes(wc2.0_bio_30s_12_Encelia, wc2.0_bio_30s_19_Encelia), col = "gray", pch = 16) +
  geom_point(data = pts4 %>% filter(flowers == "B"), aes(wc2.0_bio_30s_12_Encelia, wc2.0_bio_30s_19_Encelia), pch = 16, col = "#482624") +
  # scale_color_manual(values = c("#F5B203", "#482624")) + 
  theme_classic()

ggplot() +
  geom_point(data = pts4 %>% filter(flowers != "N"), aes(wc2.0_bio_30s_12_Encelia, wc2.0_bio_30s_19_Encelia), col = "gray", pch = 16) +
  geom_point(data = pts4 %>% filter(flowers == "Y"), aes(wc2.0_bio_30s_12_Encelia, wc2.0_bio_30s_19_Encelia), pch = 16, col = "#F5B203") +
  # scale_color_manual(values = c("#F5B203", "#482624")) + 
  theme_classic()

ggplot() + 
  geom_boxplot(data = pts4 %>% filter(flowers != "N"), aes(flowers, wc2.0_bio_30s_12_Encelia)) +
  theme_classic()
ggplot() + 
  geom_boxplot(data = pts4 %>% filter(flowers != "N"), aes(flowers, wc2.0_bio_30s_19_Encelia)) +
  theme_classic()


cc %>% filter(flowers != "N") %>% group_by(flowers) %>% 
  summarize_at(grep("bio", names(cc)), mean) %>% View()
  
