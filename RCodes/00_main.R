### Cargar librerias a utiizar
library(jsonlite)
library(aqp)
library(ithir)
library(sp)
library(lattice)
#library(cluster)
library(sharpshootR)
library(sf)
library(gstat)
library(raster)
library(rgdal)
library(tmap)

# Cambiar directorio de trabajo
setwd( "D:/Tesis/Datos")

#raster::shapefile()

#cargar y leer los archivos a utilizar
aoi <- shapefile("areaEstudioUTM.shp")
json <- fromJSON("perfilesSislac.json", flatten = TRUE)$features

### OBJETIVO ESPECIFICO 1
# validacion de perfiles con error y su exclusion del archivo JSON
errores <- f0101_validarErrores(json)
json <- f0102_removerErrores(json,errores)
# Validacion de inconsistencias corregibles en los horizontes
inconsistencias <- f0103_validarInconsistencias(json)
# Corregir inconsistencias y armar objeto SoilProfileCOllection
spc <- f0105_crearSPC(json, inconsistencias)


### OBJETIVO ESPECIFICO 2
spc1 <- f0201_segmentacionAdaptada(spc)


### OBJETIVO ESPECIFICO 3
spc1 <- f0301_segmentacionAQP(spc)    # Segmentacion del spc para COS
aggAQP2 <- f0302_agregacionAQP(spc1)     # O spc2, para los valores ajustados
f0303_ploteoAgregacion(aggAQP2)          # 


### OBJETIVO ESPECIFICO 4
# Se carga la zona de estudio, shapefile en UTM
interpolacion <- f0401_interpolacionIDW(spc1, aoi) # O spc2, para los valores ajustados
