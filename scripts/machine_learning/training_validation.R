library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/machine_learning/")
x = read.csv("Encelia.itzel.csv")
y = read.csv("Encelia.sonal1.csv")
z = read.csv("Encelia.sonal2.csv")

xx = left_join(x %>% dplyr::select(file, i = flowers), y %>% dplyr::select(file, s = flowers))
for (i in 1:nrow(z)) {
  xx[which(xx$file == z[i, "file"]), "s"] = z[i, "flowers"]
}

xx1 = xx %>% filter(complete.cases(s)) %>%
  filter(i %in% c("B", "Y", "N"), s %in% c("B", "Y", "N"))
xx2 = xx1 %>% group_by(i, s) %>% summarize(n = n()) %>%
  mutate(same = ifelse(i == s, TRUE, FALSE))
xx2 %>% group_by(same) %>% summarize(n = sum(n))
xx3 = xx1 %>% filter(i == s)
write.csv(xx3, "Encelia.training_validation.csv")