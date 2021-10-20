# Cargar paquetes

library(vegan)
library(tidyverse)
library(sf)
source('biodata/funciones.R')

load('biodata/Chrysobalanaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')

# Analisis 1

#Lista de especies: 4
sort(colnames(mc_chrys))

#Abundancia de toda la comunidad
sum(colSums(mc_chrys))

#Abundancia por cuadrante
sort(rowSums(mc_chrys))

#Abundancia por especie
sort(colSums(mc_chrys))

#Tabla de abundancias
abun_sp <- censo_chrys %>%
  group_by(Latin) %>% 
  count() %>% 
  arrange(desc(n))
abun_sp

