rm(list = ls())

library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/")

y = read.csv("machine_learning/accuracy_set.csv")

a = list.files("machine_learning/categorization/", full.names = T)
lr = gsub(".csv", "", gsub(".*lr_", "", a))
res = vector("list", length(a))
for (i in 1:length(a)) {
  x = read.csv(a[i])
  xx = x %>% left_join(y)
  xx[xx$prediction == 0, "prediction"] = "B"
  xx[xx$prediction == 1, "prediction"] = "N"
  xx[xx$prediction == 2, "prediction"] = "Y"
  xx = xx %>% mutate(same = ifelse(flowers == prediction, TRUE, FALSE),
                     lr = lr[i])
  res[[i]] = xx
}

res2 = do.call("rbind", res)

# get accuracy rate including no flowers
res2 %>% group_by(lr) %>% summarize(accuracy = sum(same) / n())

# get accuracy rate for flowers only
res2 %>% filter(prediction %in% c("B", "Y")) %>% group_by(lr) %>% 
  summarize(accuracy = sum(same) / n())

# for that accuracy rate, get whether saying black when yellow or no flower when flower
res2 %>% filter(prediction %in% c("B", "Y")) %>% 
  filter(same == FALSE, lr == "0.01") %>% dplyr::select(file, prediction, flowers, X0, X1, X2)
  
