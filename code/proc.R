library(vegan)
library(tidyverse)

source("utils.R")

#data_full = load_data()
data_full = read_csv("../data/new_data.csv")
data_def_nona = data_full %>%
  select(var_def) %>% # defined in utils.R
  na.omit()
data_all_nona = data_full %>%
  select(var_all) %>% # defined in utils.R
  na.omit()

draw_climate_pca(data_full)

# to analysis with site averages, uncomment the following lines
#
data_def_avg = data_full %>%
  select(all_of(var_def)) %>%
  group_by(site) %>%
  summarise_at(var_defense, mean, na.rm=TRUE) %>%
  left_join(select(data_full, site, origin), by="site") %>%
  distinct() %>%
  na.omit()

lst = clust_and_pca(data_def_avg)

lst = clust_and_pca(data_def_nona)
draw_defense_pca(lst$st, lst$sp, lst$gr, lst$pca_defense)

aov.tab = permutest_cluster(data_def_nona, lst)
write.csv(aov.tab, "../results/permtest_cluster.csv")

#################################################################

data = centroids_and_area(data_def_nona, data_all_nona)
data$PC1 = lst$st[,1]
data$PC2 = lst$st[,2]
data_site = data %>% 
  group_by(site) %>% 
  summarise_at(c(var_defense, var_clim, "area", "cent.x", "cent.y", "herbivory", "PC1", "PC2"), mean) %>%
  left_join(select(data, site, origin), by = "site")

data_site = inner_join(data_site, select(data_full, site, latitude), by="site")

#multiple_mantel(data, "all")
#multiple_mantel(data %>% filter(origin=="native"), "native")
#multiple_mantel(data %>% filter(origin=="introduced"), "introduced")

rda_and_summ(data, var_defense, "../results/quadrat_traits.csv")
rda_and_summ(data, c("cent.x", "cent.y"), "../results/quadrat_cent.csv")
rda_and_summ(data, c("area"), "../results/quadrat_area.csv")

var = var_defense
rvar = data_site %>% select(var) %>% scale() %>% as.data.frame()
evar = cbind(
  data_site %>% select(var_clim) %>% scale(),
  data_site %>% select(var_herb) %>% scale(),
  data_site %>% select(origin, latitude)
  
) %>% as.data.frame()

ret = rda(
  rvar ~
    origin*latitude,
  data = evar
)
perm <- how(nperm = 999)
#setBlocks(perm) <- with(data, site)
aov.rda = anova(ret, by="term", permutations = perm)

write.csv(aov.rda, fname)
aov.rda

rda_and_summ(data_site, var_defense, "../results/site_traits.csv")
rda_and_summ(data_site, c("cent.x", "cent.y"), "../results/site_cent.csv")
rda_and_summ(data_site, c("area"), "../results/site_area.csv")
