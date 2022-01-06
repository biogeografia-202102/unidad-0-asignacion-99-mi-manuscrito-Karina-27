#' Script Reproducible
#' 
#' 
#' Cargar paquetes
#'
library(vegan)
library(tidyverse)
library(sf)
source('biodata/funciones.R')
#'
#' Cargar datos
load('biodata/Chrysobalanaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
bci_env_grid %>% tibble
censo_chrys %>% tibble
#'
#' Analisis Exploratorio
#' Lista de especies: 4
sort(colnames(mc_chrys))
#'
#' Abundancia de toda la comunidad
sum(colSums(mc_chrys))
#'
#' Abundancia por cuadrante
sort(rowSums(mc_chrys))
#'
#' Abundancia por especie
sort(colSums(mc_chrys))
#'
#' Tabla de abundancias
abun_sp <- censo_chrys %>%
  group_by(Latin) %>% 
  count() %>% 
  arrange(desc(n))
abun_sp
#'
#' Analisis de Agrupamiento
#' Cargar otros paquetes
library(cluster)
library(gclus)
library(pvclust)
library(mapview)
library(RColorBrewer)
library(broom)
library(indicspecies)
#'
#' Crear abreviaturas de la familia
mi_fam <- mc_chrys
(colnames(mi_fam) <- make.cepnames(colnames(mi_fam)))
(df_equivalencias <- data.frame(
  nombre_original = colnames(mc_chrys),
  colnames(mi_fam)))
#'
grupos_ward_k2 <- readRDS('grupos_ward_k2.RDS')
table(grupos_ward_k2)
grupos_compl_k2 <- readRDS('grupos_compl_k2.RDS')
table(grupos_compl_k2)
#'
#'#' Cargar paletas de colores
rojo <- colorRampPalette(brewer.pal(8, "Reds"))
rojo_inv <- colorRampPalette(rev(brewer.pal(8, "Reds")))
colores_grupos <- brewer.pal(8, "Accent")
#'
#' Transformar la matriz de comunidad
#'
mi_fam_norm <- decostand(mi_fam, "normalize")
mi_fam_norm_d <- vegdist(mi_fam_norm, "euc")
mi_fam_norm_d %>% tidy
attr(mi_fam_norm_d, "labels") <- rownames(mi_fam)
#' 
#' Agrupamiento por Enlace completo
(cl_complete <- hclust(mi_fam_norm_d, method = 'complete'))
plot(cl_complete, labels = rownames(mi_fam), hang = -1,
     main = "Sitios de BCI según composición de especies de Chrysobalanaceae\nEnlace completo a partir de matriz de distancia de cuerdas",
     xlab = 'Sitios', ylab = 'Altura')
#' Agrupamiento por medio de Ward
(cl_ward <- hclust(mi_fam_norm_d, method = 'ward.D2'))
plot(cl_ward, labels = rownames(mi_fam), hang = -1,
     main = "Sitios de BCI según composición de especies de Chrysobalanaceae\nMétodo de Ward a partir de matriz de distancia de cuerdas",
     xlab = 'Sitios', ylab = 'Altura')
#' 
#' Medicion de la correlacion cofenetica
#' 
lista_cl <- list(
  cl_complete = hclust(mi_fam_norm_d, method = 'complete'),
  cl_ward = hclust(mi_fam_norm_d, method = 'ward.D2')
)
par(mfrow = c(1,2))
invisible(map(names(lista_cl), function(x) plot(lista_cl[[x]], main = x, hang = -1, cex.lab = 1)))
par(mfrow = c(1,1))
#'
map_df(lista_cl, function(x) {
  coph_d <- cophenetic(x)
  corr <- cor(mi_fam_norm_d, coph_d)
  return(corr)
})
#'
#' Calculo de la anchura de siluetas
#'
anch_sil_ward <- calcular_anchuras_siluetas(
  mc_orig = mi_fam, 
  distancias = mi_fam_norm_d, 
  cluster = lista_cl$cl_ward)
anch_sil_ward
w_dend_reord <- reorder.hclust(lista_cl$cl_ward, mi_fam_norm_d)
plot(w_dend_reord, hang = -1)
rect.hclust(
  tree = w_dend_reord,
  k = anch_sil_ward$n_grupos_optimo)
heatmap(
  as.matrix(mi_fam_norm_d),
  Rowv = as.dendrogram(w_dend_reord),
  symm = TRUE,
  margin = c(3, 3),
  col = rev(cm.colors(4))
)
#' Bootstrap
cl_pvclust_ward <-
  pvclust(t(mi_fam_norm),
          method.hclust = "ward.D2",
          method.dist = "euc",
          iseed = 191, # Resultado reproducible
          parallel = TRUE)
plot(cl_pvclust_ward, hang = -1)
lines(cl_pvclust_ward)
pvrect(cl_pvclust_ward, alpha = 0.91, border = 4)
saveRDS(grupos_ward_k2, 'grupos_ward_k2.RDS')
#' Para complete
anch_sil_compl <- calcular_anchuras_siluetas(
  mc_orig = mi_fam, 
  distancias = mi_fam_norm_d, 
  cluster = lista_cl$cl_complete)
anch_sil_compl
(grupos_compl_k2 <- as.factor(cutree(lista_cl$cl_complete, k = 2)))
table(grupos_compl_k2)
#'
w_dend_reord <- reorder.hclust(lista_cl$cl_complete, mi_fam_norm_d)
plot(w_dend_reord, hang = -1)
rect.hclust(
  tree = w_dend_reord,
  k = anch_sil_compl$n_grupos_optimo)
#'
heatmap(
  as.matrix(mi_fam_norm_d),
  Rowv = as.dendrogram(w_dend_reord),
  symm = TRUE,
  margin = c(3, 3),
  col = rev(cm.colors(4))
)
#' Bootstrap
cl_pvclust_compl <-
  pvclust(t(mi_fam_norm),
          method.hclust = "complete",
          method.dist = "euc",
          iseed = 27, # Resultado reproducible
          parallel = TRUE)
plot(cl_pvclust_compl, hang = -1)
lines(cl_pvclust_compl)
pvrect(cl_pvclust_compl, alpha = 0.91, border = 4)
(grupos_compl_k2 <- as.factor(cutree(lista_cl$cl_complete, k = 2)))
table(grupos_compl_k2)
saveRDS(grupos_compl_k2, 'grupos_compl_k2.RDS')
#'
#' Relación de Variables ambientales con grupos de Ward
(m_amb_ward_k2 <- bci_env_grid %>%
    select_if(is.numeric) %>% select(-id) %>% 
    mutate(grupos_ward_k2) %>%
    st_drop_geometry() %>% 
    pivot_longer(-grupos_ward_k2, names_to = "variable", values_to = "valor"))
#' 
m_amb_ward_k2 %>%
  group_by(variable) %>%
  summarise(
    p_valor_t = t.test(valor ~ grupos_ward_k2)$p.value,
    p_valor_w = wilcox.test(valor ~ grupos_ward_k2, exact = F)$p.value) %>%
  arrange(p_valor_t) %>%
  print(n=Inf)
#' 
m_amb_ward_k2 %>% 
  group_by(variable) %>% 
  ggplot() + aes(x = grupos_ward_k2, y = valor, fill = grupos_ward_k2) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = 'Accent') +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~ variable, scales = 'free_y')
#' 
#' Mapa de grupos por medio de Ward
mapa_ward_k2 <- mapView(
  bci_env_grid %>% mutate(grupos_ward_k2),
  layer.name = 'Grupos (2) Ward',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = colores_grupos[1:2],
  zcol = 'grupos_ward_k2') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 16)
mapa_ward_k2
#' 
#' Especies indicadoras para Ward_k2
iva_ward_k2 <- multipatt(
  x = mi_fam,
  cluster = grupos_ward_k2,
  func = 'IndVal.g',
  max.order = 2,
  control = how(nperm = 999))
summary(iva_ward_k2, indvalcomp = TRUE)
colSums(mi_fam)
(p_ward_adj <- p.adjust(iva_ward_k2$sign$p.value))
(iva_ward_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_ward_k2,
  func = "IndVal.g",
  nboot = 1000))
#'
phi_ward_k2 <- multipatt(
  mi_fam,
  grupos_ward_k2,
  func = "r.g",
  max.order = 2,
  control = how(nperm = 999))
summary(phi_ward_k2)
colSums(mi_fam)
(phi_ward_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_ward_k2,
  func = "r.g",
  nboot = 1000))
#'
#' Cargar otros paquetes
library(ape)
library(spdep)
library(ade4)
library(adegraphics)
library(adespatial)
library(gridExtra)
library(grid)
library(gtable)
library(magrittr)
library(plyr)
library(SpadeR)
library(iNEXT)
library(vegetarian)
source('https://raw.githubusercontent.com/maestria-geotel-master/unidad-3-asignacion-1-vecindad-autocorrelacion-espacial/master/lisaclusters.R')
#' Tecnicas de ordenacion
#' 
#' Hacer una matriz de variables ambientales solo con variables seleccionadas
env_select <- bci_env_grid %>% 
  st_drop_geometry %>%
  dplyr::select_if(is.numeric) %>%
  dplyr::select(-id) %>%
  dplyr::select(Al, Cu, Mn, N, N.min., pH, elevacion_media, heterogeneidad_ambiental, orientacion_media)
env_select %>% tibble
#' 
#' Analisis de Redundancia
mi_fam_hel <- decostand(mi_fam, method = 'hellinger')
mi_fam_hel %>% tibble
mi_fam_hel_rda_select <- rda(mi_fam_hel ~ ., env_select)
summary(mi_fam_hel_rda_select)
#'
RsquareAdj(mi_fam_hel_rda_select)$adj.r.squared
vif.cca(mi_fam_hel_rda_select)
#' 
#' Analisis de Componentes Principales
env_select_pca <- rda(env_select, scale = TRUE)
env_select_pca
summary(env_select_pca)
#'
par(mfrow = c(1, 2))
cleanplot.pca(env_select_pca, scaling = 1, mar.percent = 0.08, cex.char1 = 0.5)
cleanplot.pca(env_select_pca, scaling = 2, mar.percent = 0.04, cex.char1 = 0.5)
par(mfrow = c(1, 1))
#'
(env_agrupamiento <- hclust(dist(scale(env_select)), 'ward.D'))
(env_grupos <- cutree(env_agrupamiento, k = 2))
(mi_cluster <- factor(env_grupos))
(mi_cluster_l <- levels(mi_cluster))
(mi_cluster_l_seq <- 1:length(mi_cluster_l))
(puntuaciones <- scores(env_select_pca, display = 'wa', scaling = 1))
par(mfrow = c(1, 2))
grafico_base <- plot(
  env_select_pca,
  display = "wa",
  scaling = 1,
  type = "n",
  main = "PCA según grupos por variables seleccionadas"
)
abline(v = 0, lty = "dotted")
abline(h = 0, lty = "dotted")
for (i in mi_cluster_l_seq) {
  points(puntuaciones[mi_cluster == i, ],
         pch = (14 + i),
         cex = 2,
         col = i + 1)
}
text(puntuaciones, row.names(env_select), cex = 1, pos = 3)
legend(
  "topleft", # Otras alternativas: "bottomleft", "bottomright" y "topleft"
  paste("Grupo", c(mi_cluster_l_seq)),
  pch = 14 + c(mi_cluster_l_seq),
  col = 1 + c(mi_cluster_l_seq),
  pt.cex = 2
)
#'
(mi_cluster_anterior <- grupos_compl_k2)
(mi_cluster_anterior_l <- levels(mi_cluster_anterior))
(mi_cluster_anterior_l_seq <- 1:length(mi_cluster_anterior_l))
grafico_base <- plot(
  env_select_pca,
  display = "wa",
  scaling = 1,
  type = "n",
  main = "PCA según grupos por enlace completo"
)
abline(v = 0, lty = "dotted")
abline(h = 0, lty = "dotted")
for (i in mi_cluster_anterior_l_seq) {
  points(puntuaciones[mi_cluster_anterior == i, ],
         pch = (14 + i),
         cex = 2,
         col = i + 1)
}
text(puntuaciones, row.names(env_select), cex = 1, pos = 3)
legend(
  "topleft", # Otras alternativas: "bottomleft", "bottomright" y "topleft"
  paste("Grupo", c(mi_cluster_anterior_l_seq)),
  pch = 14 + c(mi_cluster_anterior_l_seq),
  col = 1 + c(mi_cluster_anterior_l_seq),
  pt.cex = 2
)
par(mfrow = c(1, 1))
#'
#' PCA a datos de comunidad
mi_fam_hel <- decostand(mi_fam, method = 'hellinger')
mi_fam_hel %>% tibble
mi_fam_hel_pca <- rda(mi_fam_hel)
summary(mi_fam_hel_pca)
#'
biplot(
  mi_fam_hel_pca,
  main = "PCA escalamiento 2 con ajuste a variables seleccionadas")
(mi_fam_hel_pca_envfit <- envfit(mi_fam_hel_pca, env_select, scaling = 2))
plot(mi_fam_hel_pca_envfit, p.max = 0.05 , col = 3)
#'
#' Analisis de Correspondencia
mi_fam_ca <- cca(mi_fam)
summary(mi_fam_ca)
summary(mi_fam_ca, scaling = 1)
par(mfrow = c(1, 2))
plot(mi_fam_ca,
     scaling = 1,
     main = "Análisis de correspondencia, escalamiento 1"
)
plot(mi_fam_ca,
     scaling = 2, # Por defecto scaling=2, lo escribo sólo para fines didáticos
     main = "Análisis de correspondencia, escalamiento 2")
par(mfrow = c(1, 1))
#'
#' Analisis de Correspondencia Canonica
mi_fam_cca_select <- cca(mi_fam ~ ., env_select)
summary(mi_fam_cca_select)
RsquareAdj(mi_fam_cca_select)
#' 
plot(mi_fam_cca_select,
     scaling = 2,
     display = c("sp", "lc", "cn"),
     main = "Triplot de CCA especies según variables seleccionadas, escalamiento 2"
)
#' 
#' Analisis de Diversidad
#' 
(indices <- alpha_div(mi_fam))
indices_sitios <- indices[-c(6,7,34,38,39),] #Para sitios con mas de una especie
pairs(indices_sitios,
      lower.panel = panel.smooth,
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      main = "Pearson Correlation Matrix")
#' 
indices_env <- bind_cols(
  indices_sitios,
  bci_env_grid[-c(6,7,34,38,39),] %>%
    select_if(is.numeric) %>%
    st_drop_geometry %>%
    select(-id) %>% 
    select(Mg, Zn, N.min., N, pH, orientacion_media, geomorf_interfluvio_pct))
indices_env %>% tibble
ezCorM(indices_env, r_size_lims = c(3,5), label_size = 4)
#'
#' Modelos de abundancia de especies
riqueza <- specnumber(mi_fam)
abundancia <- rowSums(mi_fam)
(rango_abun <- range(abundancia))
riqueza_menor_abun <- rarefy(mi_fam, sample = rango_abun[1])
sort(riqueza)
sort(round(riqueza_menor_abun))
rarecurve(
  mi_fam,
  step = 1,
  sample = rango_abun[1],
  xlab = "Número de individuos (tamaño de muestra)",
  ylab = "Especies",
  label = TRUE,
  col = "blue"
)
#' Estimacion de la riqueza para toda la comunidad y para los grupos
mi_fam_combinada <- colSums(mi_fam)
mi_fam_combinada %>% sort
mi_fam_combinada_chao <- estimacion_riqueza_chao(
  mc = mi_fam_combinada,
  n_raras = 22)
mi_fam_combinada_chao$asintoticos_estimacion
mi_fam_combinada_chao$no_asintoticos_rarefaccion_extrapolacion
mi_fam_combinada_chao$no_asintoticos_rarefaccion_extrapolacion_grafico
#'
mi_fam_k2 <- mi_fam %>%
  mutate(g=grupos_ward_k2) %>%
  group_by(g) %>%
  summarise_all(sum) %>%
  select(-g) %>% 
  data.frame
mi_fam_k2 %>% rowSums %>% sort
mi_fam_k2_chao <- estimacion_riqueza_chao(
  mc = mi_fam_k2,
  n_raras = 22)
mi_fam_k2_chao$asintoticos_estimacion
mi_fam_k2_chao$no_asintoticos_rarefaccion_extrapolacion
mi_fam_k2_chao$no_asintoticos_rarefaccion_extrapolacion_grafico
#' 
#' Diversidad beta
determinar_contrib_local_y_especie(
  mc = mi_fam,
  alpha = 0.05,
  nperm = 9999,
  metodo = 'hellinger')
#' 
#' Analisis de Ecologia Espacial
#'
#' Transformar matriz a Hellinger y generar vecindad
mi_fam_hel <- decostand (mi_fam, "hellinger")
bci_env_grid_sp <- bci_env_grid %>% as_Spatial
centroides <- bci_env_grid %>% st_centroid
bci_xy <- centroides %>% st_coordinates %>% as.data.frame
(vecindad <- bci_env_grid_sp %>% poly2nb)
(pesos_b <- nb2listw(vecindad, style = 'B'))
plot(bci_env_grid_sp)
plot(vecindad, coords = bci_xy, add=T, col = 'red')
#'
#' Autocorrelacion mediante correlograma para variables ambientales
bci_env_grid_num <- bci_env_grid %>%
  st_drop_geometry %>% 
  select_if(is.numeric) %>% 
  select(B, Ca, Mg, Zn, N, pH, N.min.)
suppressWarnings(auto_amb <- calcular_autocorrelacion(
  df_fuente = bci_env_grid_num,
  orden = 9,
  obj_vecindad = vecindad))
print(auto_amb, digits = 2, p.adj.method = 'holm')
dim_panel <- rev(n2mfrow(ncol(bci_env_grid_num)))
par(mfrow = dim_panel)
suppressWarnings(invisible(lapply(auto_amb, function(x) plot(x, main = x$var))))
#'
#' Correlograma de Mantel sin tendencias espaciales
mi_fam_sin_tendencia <- resid(
  lm(as.matrix(mi_fam_hel) ~ .,
     data = bci_xy))
mi_fam_sin_tendencia_d <- dist(mi_fam_sin_tendencia)
(mi_fam_correlograma <- mantel.correlog(
  mi_fam_sin_tendencia_d,
  XY = bci_xy,
  nperm = 999))
plot(mi_fam_correlograma)
#'
#' Calculo de I de Moran y mapas Lisa para especies sin tendencia y algunas variables ambientales
#'
(autocor_global_residuos <- sapply(
  dimnames(mi_fam_sin_tendencia)[[2]],
  function(x)
    moran.mc(
      x = mi_fam_sin_tendencia[,x],
      listw = pesos_b,
      zero.policy = T,
      nsim = 9999),
  simplify = F))
#'
bci_env_grid_num_sf <- bci_env_grid %>%
  select_if(is.numeric) %>% 
  select(-id, -UTM.EW, -UTM.NS)
bci_env_grid_num_sf %>% tibble
lisamaps_amb <- sapply(grep('geometry', names(bci_env_grid_num_sf), invert = T, value = T),
                       function(x) {
                         m <- lisamap(objesp = bci_env_grid_num_sf[x],
                                      var = x,
                                      pesos = pesos_b,
                                      tituloleyenda = 'Significancia ("x-y", léase como "x" rodeado de "y")',
                                      leyenda = F,
                                      anchuratitulo = 50,
                                      tamanotitulo = 10,
                                      fuentedatos = '\nhttp://ctfs.si.edu/webatlas/datasets/bci/',
                                      titulomapa = paste0('Clusters LISA de "', x, '"'))
                         return(m$grafico)
                       }, simplify = F
)
lisamaps_amb$leyenda <- gtable_filter(ggplot_gtable(ggplot_build(lisamaps_amb[[1]] + theme(legend.position="bottom"))), "guide-box")
grid.arrange(do.call('arrangeGrob', c(lisamaps_amb[1:12], nrow = 3)), lisamaps_amb$leyenda, heights=c(1.1, 0.1), nrow = 2)
grid.arrange(do.call('arrangeGrob', c(lisamaps_amb[13:22], nrow = 3)), lisamaps_amb$leyenda, heights=c(1.1, 0.1), nrow = 2)
grid.arrange(do.call('arrangeGrob', c(lisamaps_amb[23:31], nrow = 3)), lisamaps_amb$leyenda, heights=c(1.1, 0.1), nrow = 2)
#'
mi_fam_sintendencia_sf <- bci_env_grid %>% select %>% bind_cols(mi_fam_sin_tendencia %>% as.data.frame)
lisamaps_mifam_sintendencia <- sapply(
  grep('geometry', names(mi_fam_sintendencia_sf), invert = T, value = T),
  function(x) {
    m <- lisamap(objesp = mi_fam_sintendencia_sf[x],
                 var = x,
                 pesos = pesos_b,
                 tituloleyenda = 'Significancia ("x-y", léase como "x" rodeado de "y")',
                 leyenda = F,
                 anchuratitulo = 50,
                 tamanotitulo = 10,
                 fuentedatos = '\nhttp://ctfs.si.edu/webatlas/datasets/bci/',
                 titulomapa = paste0('Clusters LISA de "', x, '"'))
    # dev.new();print(m$grafico)
    return(m$grafico)
  }, simplify = F
)
lisamaps_mifam_sintendencia$leyenda <- gtable_filter(ggplot_gtable(ggplot_build(lisamaps_mifam_sintendencia[[2]] + theme(legend.position="bottom"))), "guide-box")
grid.arrange(do.call('arrangeGrob', c(lisamaps_mifam_sintendencia, nrow = 3)), heights=c(1.1, 0.1))
