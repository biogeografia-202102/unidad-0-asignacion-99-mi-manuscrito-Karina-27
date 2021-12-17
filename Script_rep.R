# Cargar paquetes

library(vegan)
library(tidyverse)
library(sf)
source('biodata/funciones.R')

load('biodata/Chrysobalanaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
mi_fam <- mc_chrys
(colnames(mi_fam) <- make.cepnames(colnames(mi_fam)))
(df_equivalencias <- data.frame(
  nombre_original = colnames(mc_chrys),
  colnames(mi_fam)))
bci_env_grid %>% tibble
#'
grupos_ward_k2 <- readRDS('grupos_ward_k2.RDS')
table(grupos_ward_k2)
grupos_compl_k2 <- readRDS('grupos_compl_k2.RDS')
table(grupos_compl_k2)
#'
#' Analisis 1
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
#' Tecnicas de ordenacion
env_select <- bci_env_grid %>% 
  st_drop_geometry %>%
  dplyr::select_if(is.numeric) %>%
  dplyr::select(-id) %>%
  dplyr::select(Al, Cu, Mn, N, N.min., pH, elevacion_media, heterogeneidad_ambiental, orientacion_media)
env_select %>% tibble
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
#' Analisis de Redundancia
mi_fam_hel <- decostand(mi_fam, method = 'hellinger')
mi_fam_hel %>% tibble
mi_fam_hel_rda_select <- rda(mi_fam_hel ~ ., env_select)
summary(mi_fam_hel_rda_select)
#'
RsquareAdj(mi_fam_hel_rda_select)$adj.r.squared
vif.cca(mi_fam_hel_rda_select)
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