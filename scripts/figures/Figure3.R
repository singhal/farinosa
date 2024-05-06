rm(list = ls())
library(hzar)
library(tidyverse)
library(raster)
library(pcaMethods)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/")
x = read.csv("metadata/Encelia_Samples - GENERAL.csv")
d = read.csv("metadata/sample_colors.csv")

testing_clines <- function(cl1, cl2,  cl1_con, cl2_con, type) {
  # test for concordance and coincidence
  # by looking at AIC values
  aic1 = hzar.AIC.hzar.dataGroup(cl1)
  aic1_con = hzar.AIC.hzar.dataGroup(cl1_con)
  aic2 = hzar.AIC.hzar.dataGroup(cl2)
  aic2_con = hzar.AIC.hzar.dataGroup(cl2_con)
  print(paste("AIC for original cline model 1: ", aic1, sep=""))
  print(paste("AIC for ", type, " cline model 1: ", aic1_con, sep=""))
  aics = c(aic1, aic1_con)
  diff = exp((min(aics) - max(aics))/2)
  if (aic1_con < aic1) {
    diff = 1 - diff
  }
  
  print(paste("***The constrained model explains ", diff, " of the information.", sep=""))
  print(paste("AIC for original cline model 2: ", aic2, sep=""))
  print(paste("AIC for ", type, " cline model 2: ", aic2_con, sep=""))
  aics = c(aic2, aic2_con)
  diff = exp((min(aics) - max(aics))/2)
  
  if (aic2_con < aic2) {
    diff = 1 - diff
  }
  
  print(paste("***The constrained model explains ", diff, " of the information.", sep=""))
}

make_cline <- function(cl1) {
  clModel <- hzar.makeCline1DFreq(cl1, scaling="fixed",tails="none");
  clModel <- hzar.model.addBoxReq(clModel, 0, 1600);
  
  cl_model_fitR = hzar.first.fitRequest.old.ML(clModel, cl1, verbose=TRUE)
  
  cl_model_fitR$mcmcParam$chainLength <- 1e5;
  cl_model_fitR$mcmcParam$burnin <- 1e2;
  clModelFit2 <- hzar.doFit(cl_model_fitR)
  clModelFitData <- hzar.dataGroup.add(clModelFit2)
  clModelFitData <- hzar.dataGroup.add(clModelFitData, 
                                       hzar.chain.doSeq(hzar.next.fitRequest(clModelFit2)));
  pts = seq(0, max(cl1$frame$dist), 1)
  fzCoor <- hzar.getCredParamRed(clModelFitData)$fzCline(pts)
  predval = data.frame(distance = pts,
                       percent_black = clModelFitData$ML.cline$clineFunc(pts))
  return(list(clModelFitData, predval))
}

make_cline2 <- function(cl1, to_fix = NULL, val = NULL) {
  clModel <- hzar.makeCline1DFreq(cl1, scaling="fixed",tails="none");
  clModel <- hzar.model.addBoxReq(clModel, 0, 1600);
  
  if (!missing(to_fix)) {
    clModel[['parameterTypes']][[to_fix]]$val <- val
    attr(clModel[['parameterTypes']][[to_fix]], 'fixed') <- TRUE
  }
  
  cl_model_fitR = hzar.first.fitRequest.old.ML(clModel, cl1, verbose=TRUE)
  
  cl_model_fitR$mcmcParam$chainLength <- 1e5;
  cl_model_fitR$mcmcParam$burnin <- 1e2;
  clModelFit2 <- hzar.doFit(cl_model_fitR)
  clModelFitData <- hzar.dataGroup.add(clModelFit2)

  return(clModelFitData)
}

##############
# color cline
##############

c = read.csv("metadata/black_yellow.cleaned11Sept23.csv")
coors = c %>% dplyr::select(longitude, latitude)
#distance matrix
dist.in.km.matrix <- fields::rdist.earth(as.data.frame(coors), miles = F, R=6371)
#clustering
fit <- hclust(as.dist(dist.in.km.matrix), method = "single")
# group by 1 vs 5 km
c = c %>% mutate(cluster1 = cutree(fit, h = 1), 
                   cluster2 = cutree(fit, h = 5))
length(unique(c$cluster1))
length(unique(c$cluster2))

c5 = c %>% group_by(cluster2) %>% summarize(yellow=sum(yellow, na.rm=T), 
                                                                   black = sum(black, na.rm=T), 
                                                                   latitude = mean(latitude), longitude = mean(longitude)) %>% 
  ungroup() %>% mutate(total = yellow + black) %>%
  mutate(percent_black = black / total)

start = as.numeric(as.character(c5 %>% slice_max(latitude, with_ties = F) %>%
  dplyr::select(longitude, latitude)))
# distance in km
c5$distance = geosphere::distHaversine(c5[, c("longitude", "latitude")], start) / 1000

cl1 = hzar.doMolecularData1DPops(c5$distance, c5$percent_black, c5$total)
cl1res = make_cline(cl1)
range(unlist(lapply(hzar.getCredParamRed(cl1res[[1]])$clines, 
                    function(x) x$param.all$width)))

##############
# genetic cline
##############

a = read.table("admixture/encelia.samtools.miss0.7.MAC2.thinned.2.Q")
n = read.table("input_files/encelia.samtools.miss0.7.MAC2.thinned.fam")
a$ind = n$V2

a = a[!a$ind %in% c("FAR09", "FAR10"), ]
a$LATITUDE = x[match(a$ind, x$PLANT_ID), "LATITUDE"]
a$LONGITUDE = x[match(a$ind, x$PLANT_ID), "LONGITUDE"]

a2 = a %>% dplyr::select(LONGITUDE, LATITUDE)
#distance matrix
dist.in.km.matrix2 <- fields::rdist.earth(as.data.frame(a2), miles = F, R=6371)
#clustering
fit <- hclust(as.dist(dist.in.km.matrix2), method = "single")
# group by 1 vs 5 km
a2 = a %>% mutate(cluster1 = cutree(fit, h = 1), 
                 cluster2 = cutree(fit, h = 5))
length(unique(a2$cluster1))
length(unique(a2$cluster2))

a3 = a2 %>% group_by(cluster1) %>% 
  summarize(population2 = mean(V2), n = n(), latitude = mean(LATITUDE), longitude = mean(LONGITUDE)) %>% 
  ungroup()

a3$distance = geosphere::distHaversine(a3[, c("longitude", "latitude")], start) / 1000

cl2 = hzar.doMolecularData1DPops(a3$distance, a3$population2, a3$n)
cl2res = make_cline(cl2)
range(unlist(lapply(hzar.getCredParamRed(cl2res[[1]])$clines, 
                    function(x) x$param.all$width)))

##############
# eco cline
##############

climr = list.files("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/spatial_rasters/wc2.0_30s_bio/", pattern=".tif", full.names = T)
climr = stack(lapply(climr, raster))

pts = c[, c("longitude", "latitude")]
clim = raster::extract(climr, pts)

clim2 = clim[complete.cases(clim), ]
pts2 = pts[complete.cases(clim),]

cc = cbind(clim2, pts2)
xpca = pca(clim2, scale="uv", center = T, nPcs = 10)

summary(xpca)
sort(abs(xpca@loadings[, "PC1"]))
rev(sort(abs(xpca@loadings[, "PC2"])))

#distance matrix
dist.in.km.matrix3 <- fields::rdist.earth(pts2, miles = F, R=6371)
#clustering
fit3 <- hclust(as.dist(dist.in.km.matrix3), method = "single")
# group by 1 vs 5 km
pts3 = pts2 %>% mutate(cluster1 = cutree(fit3, h = 1), 
                       cluster2 = cutree(fit3, h = 5))
pts3 = cbind(pts3, xpca@scores[ , 1:2])

# explore raw data
colnames(clim) = gsub("wc2.0_bio_30s_", "", gsub("_Encelia", "", names(climr)))
r = c %>% mutate(per_black = black / (yellow+black)) %>%
  dplyr::select(per_black, cluster1) %>%
  cbind(clim) %>% 
  group_by(cluster1) %>% summarise_all(mean)

bc = read.csv("bioclim.csv")

order2 =  gsub("wc2.0_bio_30s_", "",
              gsub("_Encelia", "",
              names(rev(sort(abs(xpca@loadings[, "PC2"]))))))
order3 = bc[match(as.numeric(order2), bc$var), "meaning"]

res = vector("list", length = 18)
r2 = r %>% gather(bioclim, value, -cluster1, -per_black)
r2$bioclim = bc[match(as.numeric(r2$bioclim), bc$var), "meaning"]
r2$bioclim2 = factor(r2$bioclim, levels = order3)

# for (i in 1:length(order2)) {
#   var = order2[i]
#   res[[i]] = ggplot(r, aes(.data[[var]], per_black)) + 
#     geom_point() + theme_classic() +
#     ylab("% black")
# }
# plts = cowplot::plot_grid(plotlist = res, ncol = 5)
rr = ggplot(r2, aes(value, per_black)) +
  geom_point() + 
  ylab("% brown") + 
  facet_wrap(~bioclim2, ncol = 4, scales = "free") +
  theme_bw() + 
  theme(strip.text = element_text(size = 9, margin = margin(t = 0, r = 0, b = 2, l = 0)), 
        strip.background = element_blank())
cowplot::save_plot("figures/FigureSXX_bioclim.png", rr, base_width = 8, base_height = 10)

library(corrplot)
M = cor(r, use = "pairwise.complete.obs")
corrplot.mixed(M, order = 'AOE')

pts4 = pts3 %>% group_by(cluster2) %>% 
  summarize(PC1 = mean(PC1), PC2 = mean(PC2),  
            n = n(), latitude = mean(latitude), 
            longitude = mean(longitude)) %>% 
  ungroup()

pts4$distance = geosphere::distHaversine(pts4[, c("longitude", "latitude")], start) / 1000
pts4$PC2_normal = (pts4$PC2-min(pts4$PC2))/(max(pts4$PC2)-min(pts4$PC2))
pts4$PC1_normal = (pts4$PC1-min(pts4$PC1))/(max(pts4$PC1)-min(pts4$PC1))

supfig = ggplot() + geom_point(data = pts4, aes(distance, PC1_normal), col = "red") +
  geom_line(data = cl1res[[2]], aes(distance, percent_black), col = "#482624") +
  theme_classic() + xlab("distance (km)") + 
  # Custom the Y scales:
   scale_y_continuous(
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*1, name="% brown")
  ) +
  ylab("normalized PC1") +
  theme(
  axis.title.y = element_text(color = "black", size=13),
  axis.title.y.right = element_text(color = "red", size=13)
  ) 
cowplot::save_plot("figures/SupplFigX_PC1climate.png", supfig)
  

cl3 = hzar.doMolecularData1DPops(pts4$distance, pts4$PC2_normal, pts4$n)
cl3res = make_cline(cl3)
range(unlist(lapply(hzar.getCredParamRed(cl3res[[1]])$clines, 
                    function(x) x$param.all$width)))


cl1res[[2]]$type = "% brown flower"
cl2res[[2]]$type = "% gen. pop. 2"
cl3res[[2]]$type = "PC2 climate"
predval = rbind(cl1res[[2]], cl2res[[2]], cl3res[[2]])

cx = c5 %>% mutate(type = "% brown flower") %>% 
  dplyr::select(distance = distance, type = type, value = percent_black)
ax = a3 %>% mutate(type = "% gen. pop. 2") %>% 
  dplyr::select(distance = distance, type = type, value = population2)
px = pts4 %>% mutate(type = "PC2 climate") %>% 
  dplyr::select(distance = distance, type = type, value = PC2_normal)
ptsx = rbind(cx, ax, px)


plt = ggplot() +
  geom_point(data = ptsx, aes(distance, value, color = type), alpha = 0.5, shape = 16) +
  geom_line(data = predval, aes(distance, percent_black, color = type)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set2")) +
  xlab("distance (km)") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.2))
cowplot::save_plot("figures/Figure3.png", plt, base_height = 3, base_width = 4)


#################
# test for coincidence (same center)
#################

# flower & genetics
center_flower = as.numeric(cl1res[[1]]$ML.cline$param.free['center'])
center_gen = as.numeric(cl2res[[1]]$ML.cline$param.free['center'])
flower_gen_coi1 = make_cline2(cl1, "center", center_gen)
flower_gen_coi2 = make_cline2(cl2, "center", center_flower)
testing_clines(cl1res[[1]], cl2res[[1]], flower_gen_coi1, flower_gen_coi2, "coincidence")

# flower & climate
center_clim = as.numeric(cl3res[[1]]$ML.cline$param.free['center'])
flower_clim_coi1 = make_cline2(cl1, "center", center_clim)
flower_clim_coi2 = make_cline2(cl3, "center", center_flower)
testing_clines(cl1res[[1]], cl3res[[1]], flower_clim_coi1, flower_clim_coi2, "coincidence")


# #################
# # test for concordance (same width)
# #################
# 
# # flower & genetics
# width_flower = as.numeric(cl1res[[1]]$ML.cline$param.free['width'])
# width_gen = as.numeric(cl2res[[1]]$ML.cline$param.free['width'])
# flower_gen_con1 = make_cline2(cl1, "width", width_gen)
# flower_gen_con2 = make_cline2(cl2, "width", width_flower)
# testing_clines(cl1res[[1]], cl2res[[1]], flower_gen_con1, flower_gen_con2, "concordance")
# 
# # flower & climate
# width_clim = as.numeric(cl3res[[1]]$ML.cline$param.free['center'])
# flower_clim_con1 = make_cline2(cl1, "width", center_clim)
# flower_clim_con2 = make_cline2(cl3, "width", center_flower)
# testing_clines(cl1res[[1]], cl3res[[1]], flower_clim_con1, flower_clim_con2, "concordance")
