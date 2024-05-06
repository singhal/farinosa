library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)
library(dplyr)

world <- ne_countries(scale = "medium", returnclass = "sf")

# ff = list.files("~/Dropbox (Personal)/Encelia/analysis/black_yellow/observations//", full.names = T)
# ff = lapply(ff, read.csv)

d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/machine_learning/inat_observations/observations-391086.csv")

# d = rbind(ff[[1]], ff[[2]][, names(ff[[1]])])
# d = rbind(d, ff[[3]][, names(d)])
d1 = d[grep("inaturalist-open-data.s3.amazonaws", d$image_url), ]
d1 = d1[complete.cases(d1$latitude), ]

ggplot(data = world) +
  geom_sf() +
  geom_point(data = d1,
             aes(x = longitude,
                 y = latitude),
             shape = 16, alpha = 0.2) +
  coord_sf(xlim = c(-122, -107), ylim = c(14, 39)) +
  theme_minimal()

# note that there are some duplicates
for (i in 1:nrow(d1)) {
  id = d1[i, "id"]
  url = d1[i, "image_url"]
  newfile = paste0("/Users/singhal/Desktop/photos/", id, ".jpg")
  if (!file.exists(newfile)) {
    download.file(url, newfile, mode = "wb")
  }
}


