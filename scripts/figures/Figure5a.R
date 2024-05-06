rm(list = ls())

library(ape)
library(phytools)
library(ggplot2)
library(sf)
library(cowplot)
library(gridGraphics)
library("rnaturalearth")
library("rnaturalearthdata")

# read in species tree
t = read.tree("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
outs = c("Xylorhiza_tortifolia", "Enceliopsis_covillei")
t1 = drop.tip(t, outs)

x = read.csv("~/Desktop/flower_color.csv")
colors = x$color
names(colors) = x$sps
colors = as.factor(colors)

map = make.simmap(t1, colors, model="ER",
                       nsim=100)

cols = c("#F5B203", "#482624")
names(cols) = c("yellow", "brown")
plot_tree_func <- function() {
  plotTree(t1,ftype="i",fsize=0.5,offset=0.5,
           lwd=6)
  par(fg="transparent",lend=1)
  plotTree(t1,ftype="i",fsize=0.5,offset=0.5,
           lwd=4,color="white",add=TRUE)
  ## now plot our 100 stochastic map trees
  ## with 99% transparency
  for(i in 1:length(map)) plot(map[[i]],
                                     colors=sapply(cols,make.transparent,alpha=0.01),
                                     add=TRUE,lwd=4,ftype="i",fsize=0.5,offset=0.5)
  par(fg="black")
  nodelabels(pie=summary(map)$ace, piecol=rev(cols),
             cex=0.6)
}
pdf("~/Desktop/phylogen.pdf", width = 3, height = 2.5)
plot_tree_func()
dev.off()

# make map
xx = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")
# remove outs
xx1 = xx[grep("Encelia", xx$species), ]
xx1$species = gsub(" ", "_", xx1$species)
xx1[xx1$species == "Encelia_actonii", "species"] = "Encelia_actoni"
xx1[xx1$species == "Encelia_californica", "species"] = "Encelia_californica1"
xx1[xx1$species == "Encelia_frutescens", "species"] = "Encelia_frutescens_frutescens"
xx1$color = x[match(xx1$species, x$sps), "color"]

world <- ne_countries(scale = "medium", returnclass = "sf")
xx2 = xx1[complete.cases(xx1$color), ]
pts <- st_as_sf(xx2,
                coords = c("Longitude", "Latitude"),
                crs = st_crs(world))

m1 = ggplot() +
  geom_sf(data = world, fill = "gray95") +
  geom_sf(data = pts, aes(col = color), shape = 16) +
  coord_sf(xlim = c(-120, -108), ylim = c(23, 37)) +
  scale_color_manual(values = c("#482624", "#F5B203")) +
  theme_minimal() +
  theme(legend.position = "none")
save_plot("~/Desktop/map.pdf", m1, base_height = 4, base_width = 5)
