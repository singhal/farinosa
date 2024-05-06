rm(list = ls())
library(adegenet)
library(ggplot2)
library(patchwork)
library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/")
x = read.csv("metadata/Encelia_Samples - GENERAL.csv")
d = read.csv("metadata/sample_colors.csv")

##############
# PCA
##############


t = read.snp("input_files/encelia.samtools.miss0.7.MAC2.thinned.snp")
t <- t[!indNames(t) %in% c("FAR09", "FAR10")]

pcasnp = glPca(t, nf = 40)

inds = t@ind.names
inds = gsub(".*/", "", gsub("_merged.*", "", inds))

eig.perc <- pcasnp$eig/sum(pcasnp$eig)
pcal = data.frame(pcasnp$scores)
pcal$inds = inds

pcal$col = d[match(pcal$inds, d$individual), "color"]
pcal = cbind(pcal, x[match(pcal$inds, x$PLANT_ID), c("LATITUDE", "LONGITUDE")])

# plot(pcal$PC1, pcal$PC2, col= pcal$col)
# plot(pcal$LATITUDE, pcal$PC1, col= pcal$col)

pcal$location = ifelse(pcal$LATITUDE > 30.4 & pcal$LATITUDE < 31.7, "transition", "out")

pcaplt = ggplot(pcal, aes(PC1, PC2)) + 
  ggforce::geom_mark_ellipse(aes(x = PC1, y = PC2, fill = col), 
                             expand = unit(2, "mm"), size = 0.1, alpha = 0.3) +
  geom_point(aes(shape = location, fill = col), size = 2.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("#482624", "#F5B203")) +
  xlab(paste0("PC1 (", round(eig.perc[1] * 100, 1), "%)")) +
  ylab(paste0("PC1 (", round(eig.perc[2] * 100, 1), "%)"))


##############
# IBD
##############

get_dist <- function(row) {
  pt1 = ind[row[2], c("LONGITUDE", "LATITUDE")]
  pt2 = ind[row[3], c("LONGITUDE", "LATITUDE")]
  return(geosphere::distHaversine(pt1, pt2))
}

dd = read.csv("IBD/Encelia_farinosa_divergence.csv")
col1 = d[ match(dd$ind1, d$individual), "color"]
col2 = d[ match(dd$ind2, d$individual), "color"]
dd$type = "yellow-brown"
dd[which(col1 == "black" & col2 == "black"), "type"] = "brown-brown"
dd[which(col1 == "yellow" & col2 == "yellow"), "type"] = "yellow-yellow"

ind = x
rownames(ind) = ind$PLANT_ID
dd$geo_dist = apply(dd, 1, get_dist)
dd$inv_fst = dd$fst / (1 - dd$fst)

dd1 = dd %>% filter(fst_denom > 100, complete.cases(geo_dist), fst > -0.3) %>%
  mutate(log_dist = log(geo_dist))

inds = unique(c(dd1$ind1, dd1$ind2))
f = matrix(NA, nrow = length(inds), ncol = length(inds))
g = matrix(NA, nrow = length(inds), ncol = length(inds))
for (i in 1:length(inds)) {
  for (j in 1:length(inds)) {
    res = dd1 %>% filter(ind1 == inds[i] & ind2 == inds[j] | ind1 == inds[j] & ind2 == inds[i])
    if (nrow(res) > 0) {
      f[i, j] = res$inv_fst
      g[i, j] = res$log_dist
      }
  }
}

g[!is.finite(g)] = NA

vegan::mantel(as.dist(f, upper = T, diag = T), 
              as.dist(g, upper = T, diag = T),
              permutations = 999, na.rm = T)

ibdplt = ggplot(dd1, aes(geo_dist, inv_fst, fill = type, col = type)) + 
  geom_smooth(method = "lm") +
  geom_point(shape = 21, stroke = 0.2, col = "black") +
  theme_classic() + scale_x_log10() +
  xlab("geographic distance") + 
  ylab(expression(inverse ~ F[ST])) +
  xlim(200, 380523) +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.9),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent', color = 'transparent')) +
  scale_fill_manual(values = c("yellow-brown" = "gray", "brown-brown" = "#482624", "yellow-yellow" = "#F5B203")) +
  scale_color_manual(values = c("yellow-brown" = "gray", "brown-brown" = "#482624", "yellow-yellow" = "#F5B203"))


##############
# admixture
##############

a = read.table("admixture/encelia.samtools.miss0.7.MAC2.thinned.2.Q")
n = read.table("input_files/encelia.samtools.miss0.7.MAC2.thinned.fam")
a$ind = n$V2

a = a[!a$ind %in% c("FAR09", "FAR10"), ]

xx = tidyr::gather(a, "pop", "prob", -ind)

a = cbind(a, x[match(a$ind, x$PLANT_ID), c("LATITUDE", "LONGITUDE")])
a$color = d[match(a$ind, d$individual), "color"]

start = as.numeric(as.character(a %>% slice_max(LATITUDE, with_ties = F) %>%
                                  dplyr::select(LONGITUDE, LATITUDE)))
# distance in km
a$distance = geosphere::distHaversine(a[, c("LONGITUDE", "LATITUDE")], start) / 1000
a = a %>% arrange(distance) 

order = a %>% select(ind) %>% pull()
pcal[match(a$ind, rownames(pcal)), "location"]
coldf = data.frame(
  y = 1.03,
  x = seq(1, length(order)), 
  color = a %>% select(color) %>% pull(),
  location = pcal[match(a$ind, rownames(pcal)), "location"]
)


axplt = ggplot() +
  geom_col(data = xx, aes(factor(ind, levels = order), prob, fill = pop), color = "gray", size = 0.1) +
  geom_point(data = coldf, aes(x, y, color = color, shape = location)) + 
  theme_minimal() + labs(x = "",  y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  ylim(0, 1.04) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )  +
  scale_fill_manual(values = c("V2" = "#482624", "V1" = "#F5B203")) + 
  scale_color_manual(values = c("black" = "#482624", "yellow" = "#F5B203")) +
  scale_shape_manual(values = c("transition" = 17, "out" = 16))

layout <- "
#B
AB
AC
#C
"
pp = pcaplt + ibdplt + axplt + 
  plot_layout(design = layout, heights = c(1, 2, 2, 1), widths = c(2, 3)) +  
  plot_annotation(tag_levels = 'A')

cowplot::save_plot("figures/Figure2.pdf", pp, base_height = 6, base_width = 10)
cowplot::save_plot("figures/Figure2.png", pp, base_height = 6, base_width = 10)
