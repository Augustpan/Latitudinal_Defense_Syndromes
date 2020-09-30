# definitions of groups of variables:
var_smd = c("origin", "site", "quadrat", "latitude", "longitude")
var_defense = c("tannin", "lignin", "saponins", "flavonoid", "trichome", "cn", "sla", "ssl")
var_clim = c("MAT", "MAXT", "MINT", "ST", "AP", "WP", "DP", "SP")
var_clim = c("CPC1", "CPC2")
var_symb = c("infection", "shannon")
var_symb = c()
var_herb = c("herbivory")
var_all = c(var_smd, var_defense, var_clim, var_symb, var_herb)
var_def = c(var_smd, var_defense)

load_data = function() {
  col_spec = cols(
    .default = col_double(),
    origin = col_character(),
    site = col_character(),
    quadrat = col_character()
  )
  
  data_native = read_csv("../data/data_native.csv", col_types = col_spec)
  data_intro = read_csv("../data/data_intro.csv", col_types = col_spec)
  data_full = bind_rows(data_native, data_intro)
  
  data_full = data_full %>% mutate(ssl=.$stemlength/.$stembiomass)
  
  tpc = rda(data_full %>% select(c("MAT","MAXT","MINT","ST")) %>% scale())
  ppc = rda(data_full %>% select(c("AP","WP","DP","SP")) %>% scale())
  cpc = rda(data_full %>% select(c("MAT", "MAXT", "MINT", "ST", "AP", "WP", "DP", "SP")) %>% scale())
  
  data_full$TPC = tpc$CA$u[,1]
  data_full$PPC = ppc$CA$u[,1]
  data_full$CPC1 = cpc$CA$u[,1]
  data_full$CPC2 = cpc$CA$u[,2]
  
  data_full
}

draw_climate_pca = function(data_full) {
  cpc = data_full %>% 
    select(c("MAT", "MAXT", "MINT", "ST", "AP", "WP", "DP", "SP")) %>% 
    scale() %>%
    rda()
  summ.cpc = summary(cpc)
  st = summ.cpc$site %>% as.data.frame()
  sp = summ.cpc$species %>% as.data.frame()
  gr = data.frame(syndrome=data_full$site, origin=data_full$origin)
  pca_defense = cpc
  fig = ggplot() +
    geom_point(aes(st$PC1, st$PC2, shape=gr$origin), size=2) +
    geom_segment(aes(x=0, y=0, xend=sp$PC1, yend=sp$PC2), 
                 arrow=arrow(angle=22.5, length=unit(0.2,"cm"), type="closed"),
                 linetype=1, size=0.5, colour = "gray") +
    geom_text(aes(sp$PC1, sp$PC2, label=row.names(sp)), hjust=0.2, vjust=1.5) + 
    labs(x=paste("PC1 (", format(100 * summary(pca_defense)$cont[[1]][2,1], digits=4), "%)", sep=""),
         y=paste("PC2 (", format(100 * summary(pca_defense)$cont[[1]][2,2], digits=4), "%)", sep="")) +
    geom_hline(yintercept=0,linetype=2,size=0.2) + 
    geom_vline(xintercept=0,linetype=2,size=0.2)+
    guides(shape=guide_legend(title="Origin"),color=guide_legend(title="Origin")) + 
    scale_shape_manual(values = c(16, 1)) +
    theme_bw()
  fig
}

draw_defense_pca = function(st, sp, gr, pca_defense) {
  ggplot() +
    geom_point(aes(st$PC1, st$PC2, color=gr$syndrome, shape=gr$origin)) +
    stat_ellipse(aes(st$PC1, st$PC2, color=gr$syndrome, group=gr$syndrome)) +
    geom_segment(aes(x=0, y=0, xend=sp$PC1, yend=sp$PC2), 
                 arrow=arrow(angle=22.5, length=unit(0.2,"cm"), type="closed"),
                 linetype=1, size=0.5, colour = "red") +
    geom_text(aes(sp$PC1, sp$PC2, label=row.names(sp)), hjust=0.2, vjust=1.5) + 
    labs(x=paste("PC1 (", format(100 * summary(pca_defense)$cont[[1]][2,1], digits=4), "%)", sep=""),
         y=paste("PC2 (", format(100 * summary(pca_defense)$cont[[1]][2,2], digits=4), "%)", sep="")) +
    geom_hline(yintercept=0,linetype=2,size=0.2) + 
    geom_vline(xintercept=0,linetype=2,size=0.2)+
    guides(shape=guide_legend(title="Origin"),color=guide_legend(title="Syndrome")) + 
    theme_bw()
}

clust_and_pca = function(data) {
  hcl = data %>%
    select(var_defense) %>%
    scale() %>%
    dist() %>%
    hclust("ward.D")
  
  gr = as.factor(cutree(hcl, 4))
  
  pca_defense = data %>%
    select(var_defense) %>%
    scale() %>% 
    rda()
  
  st = summary(pca_defense)$sites %>% as.data.frame()
  sp = summary(pca_defense)$species %>% as.data.frame()
  gr = data.frame(syndrome=gr, origin=data$origin)
  
  list(st=st, sp=sp, gr=gr, pca_defense=pca_defense)
}

permutest_cluster = function(data_def_nona, lst) {
  x = data_def_nona %>% 
    mutate(syndrome=lst$gr$syndrome)
  
  x_s12 = filter(x, syndrome %in% c(1, 2))
  x_s13 = filter(x, syndrome %in% c(1, 3))
  x_s14 = filter(x, syndrome %in% c(1, 4))
  x_s23 = filter(x, syndrome %in% c(2, 3))
  x_s24 = filter(x, syndrome %in% c(2, 4))
  x_s34 = filter(x, syndrome %in% c(3, 4))
  
  # multi-groups
  multi = adonis(select(x, var_defense) ~ syndrome, data=x)
  
  # parse-wise
  pw1 = adonis(select(x_s12, var_defense) ~ syndrome, data=x_s12)
  pw2 = adonis(select(x_s13, var_defense) ~ syndrome, data=x_s13)
  pw3 = adonis(select(x_s14, var_defense) ~ syndrome, data=x_s14)
  pw4 = adonis(select(x_s23, var_defense) ~ syndrome, data=x_s23)
  pw5 = adonis(select(x_s24, var_defense) ~ syndrome, data=x_s24)
  pw6 = adonis(select(x_s34, var_defense) ~ syndrome, data=x_s34)
  
  rbind(multi$aov.tab, pw1$aov.tab, pw2$aov.tab, pw3$aov.tab, pw4$aov.tab, pw5$aov.tab, pw6$aov.tab)
}

centroids_and_area = function(data_def_nona, data_all_nona) {
  x = data_def_nona %>% 
    mutate(syndrome=lst$gr$syndrome) %>%
    select(var_defense) %>%
    scale()
  center = sweep(x, 2, apply(x, 2, min),'-')
  R = apply(x, 2, max) - apply(x, 2, min)
  x_star = sweep(center, 2, R, "/")
  
  x_star = as_tibble(x_star)[,c(3,6,2,7,8,5,1,4)]
  
  theta = (c(1,2,3,4,5,6,7,8)-1) * 2*pi/8
  xm = sweep(x_star, 2, cos(theta), "*")
  ym = sweep(x_star, 2, sin(theta), "*")
  
  centroids = t(apply(cbind(xm, ym), 1, function(x) pracma::poly_center(x[1:8], x[9:16])))
  area = apply(cbind(xm, ym), 1, function(x) pracma::polyarea(x[1:8], x[9:16]))
  
  df = tibble(
    area = area,
    cent.x = centroids[,1],
    cent.y = centroids[,2],
    syndrome = lst$gr$syndrome,
    quadrat = data_def_nona$quadrat
  )
  
  data = left_join(data_all_nona, df, by="quadrat")
  data
}

# Mantel test
multiple_mantel = function(data, prefix) {
  dmat_defense = data %>% select(var_defense) %>% scale() %>% dist()
  dmat_cent = data %>% select(cent.x, cent.y) %>% scale() %>% dist()
  dmat_area = data %>% select(area) %>% scale() %>% dist()
  dmat_clim = data %>% select(var_clim) %>% scale() %>% dist()
  dmat_symb = data %>% select(var_symb) %>% scale()%>% dist()
  dmat_herb = data %>% select(var_herb) %>% scale() %>% dist()
  
  mm_d = multi.mantel(dmat_defense, list(clim=dmat_clim, symb=dmat_symb, herb=dmat_herb))
  mm_c = multi.mantel(dmat_cent, list(clim=dmat_clim, symb=dmat_symb, herb=dmat_herb))
  mm_a = multi.mantel(dmat_area, list(clim=dmat_clim, symb=dmat_symb, herb=dmat_herb))
  
  df1 = cbind(coeff=mm_d$coefficients, t=mm_d$tstatistic, p=mm_d$probt, rsq=mm_d$r.squared)[2:4,]
  rownames(df1) = c("CLIM", "SYMB", "HERB")
  df2 = cbind(coeff=mm_a$coefficients, t=mm_a$tstatistic, p=mm_a$probt, rsq=mm_a$r.squared)[2:4,]
  rownames(df2) = c("CLIM", "SYMB", "HERB")
  df3 = cbind(coeff=mm_c$coefficients, t=mm_c$tstatistic, p=mm_c$probt, rsq=mm_c$r.squared)[2:4,]
  rownames(df3) = c("CLIM", "SYMB", "HERB")
  write.csv(rbind(df1,df2,df3), paste0("../",prefix,"_mantel.csv"))
}

rda_and_summ = function(data, var, fname) {
  rvar = data %>% select(var) %>% scale() %>% as.data.frame()
  evar = cbind(
    data %>% select(var_clim) %>% scale(),
    data %>% select(var_herb) %>% scale(),
    data %>% select(origin)
  ) %>% as.data.frame()
  
  ret = rda(
    rvar ~
      CPC1 * CPC2 * herbivory * origin,
    data = evar
  )
  perm <- how(nperm = 999)
  #setBlocks(perm) <- with(data, site)
  aov.rda = anova(ret, by="term", permutations = perm)
  
  write.csv(aov.rda, fname)
  aov.rda
}

vp_and_plot_EOL = function(data, var, fname) {
  rvar = data %>% select(var) %>% scale() %>% as.data.frame()
  env = cbind(
    data %>% select(var_clim) %>% scale(),
    data %>% select(var_symb) %>% scale(),
    data %>% select(var_herb) %>% scale()
  )
  org = data %>% select(origin)
  lat = data %>% select(latitude)
  
  vp = varpart(rvar, env, lat, org)
  plot(vp, Xnames=c("ENV","LAT","ORG"))
  title(fname)
}

vp_and_plot_ABO = function(data, var, fname) {
  rvar = data %>% select(var) %>% scale() %>% as.data.frame()
  bio = cbind(
    data %>% select(var_symb) %>% scale(),
    data %>% select(var_herb) %>% scale()
  )
  abio = data %>% select(var_clim) %>% scale()
  org = data %>% select(origin)
  #lat = data %>% select(latitude)
  
  vp = varpart(rvar, bio, abio, org)
  plot(vp, Xnames=c("BIO","ABIO","ORG"))
  title(fname)
}