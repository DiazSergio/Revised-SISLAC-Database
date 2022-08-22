# We will need some packages for (spatial) data processing
library(tidyverse) # wrangling tabular data and plotting
library(sf) # processing spatial vector data - the easy way
library(sp) # processing spatial vector data - the way gstat needs it
library(raster) # processing spatial raster data. !!!overwrites dplyr::select!!!

# Packages for geostatistics
library(gstat)   # The most popular R-Package for Kriging (imho)
library(automap) # Automatize some (or all) parts of the gstat-workflow 

# Finally, some packages to make pretty plots
library(patchwork)
library(viridis)
library(tmap)

setwd("D:/Tesis/Datos/Kriging")

aoi <- shapefile("areaEstudioUTM.shp")
sitios <- readr::read_csv("Kriging1.csv") #%>% dplyr::select(-licence)

#Histograma, si no hay distribución normal se aplica transformación logarítmica.
hist(sitios$p0a5, breaks = 16)
#sitios$p0a5 <- sitios$p0a5 
summary(sitios$p0a5)
sitios$p0a5 <- (sitios$p0a5) + 1
sitios$p0a5 <- log10(sitios$p0a5)
hist(sitios$p0a5, breaks = 16)

#Calcule cuántos pares de puntos hay en el dataset meuse.
n <- length(sitios$p0a5)
n * (n - 1)/2

sitiosSF<- st_as_sf(sitios, coords = c("X", "Y"), crs = 4326)
sitiosMagna <- st_transform(sitiosSF, proj4string(aoi)) #En este caso CRS 3116, Colombia Bogota
sitiosCOS <- as(sitiosMagna, 'Spatial')
sitiosCOS@bbox <- aoi@bbox

#v_emp_OK <- gstat::variogram(p0a5~1, as(wellobs_sf, "Spatial"))
v_emp_OK <- gstat::variogram(p0a5~1, sitiosCOS)
#plot(v_emp_OK)

# automap's autofitVariogram actually produces more info than we need.
# I will only keep the var_model part.
v_mod_OK <- automap::autofitVariogram(p0a5~1, sitiosCOS)$var_model
plot(automap::autofitVariogram(p0a5~1, sitiosCOS))

#GRILLA
# Se debe calcular con base en el bbox de aio el número de celdas para una escala aprox. de 250m
nCels <- as.integer(((aoi@bbox[1,2]-aoi@bbox[1,1])*(aoi@bbox[2,2] - aoi@bbox[2,1]))/62500)
grd <- as.data.frame(spsample(sitiosCOS, "regular", n=nCels))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(sitiosCOS)

#Kriging
OK <- krige(p0a5~1, sitiosCOS, grd, model = v_mod_OK)
OK$var1.pred<- 10^(OK$var1.pred)#volver a valores originales
#OK$var1.pred<- OK$var1.pred
str(OK)

rNew <- mask(raster(OK), aoi)

titleMap <- paste("SOC (%) Depth: 0-5 cm")
breaks = c(0, 1.20, 2.40,20.00)

mapaOK <- tm_shape(rNew) + tm_raster(palette = "YlOrBr", title= titleMap, breaks= breaks,
                                  labels = c("< 1.2", "1.2 - 2.4", "> 2.4")) +
  tm_legend(outside = FALSE, text.size = .8, legend.title.fontface = "bold") +  
  tm_layout(frame = TRUE, inner.margins = .15) + 
  tm_scale_bar(position=c("left", "bottom"), breaks = c(0, 10, 20, 30), text.size = 1) + 
  tm_grid(labels.inside.frame = FALSE, n.x = 3, n.y = 4, col='gray71') +
  tm_shape(aoi) +  tm_borders()


#Cross validation
IDW.out <- vector(length = length(sitiosCOS))
for (k in 1:length(sitiosCOS)) {
  IDW.out[k] <- krige(p0a5~1, sitiosCOS[-k,], sitiosCOS[k,], model = v_mod_OK)$var1.pred
}

IDW.out <- 10^(IDW.out)

# RMSE
valorRMSE <- sqrt(sum((IDW.out - sitiosCOS$p0a5)^2) / length(sitiosCOS))
valorRMSE <- 10^(valorRMSE)
#idwRMSE <- c(idwRMSE, valorRMSE) # Tercer objeto de la lista
#R2
rss <- sum((sitiosCOS$p0a5 - IDW.out)^2)
rss <- 10^(rss)
tss <- sum((sitiosCOS$p0a5 - mean(sitiosCOS$p0a5))^2)
tss <- 10^(tss)
valorR2 <-  1 - (rss/tss)
valorR2 <- 10^(valorR2)

titleGraph <- paste("Cross validation, 15-30 cm " ,
                    "\nRMSE:", round(valorRMSE,2),
                    "- R2:",round(valorR2,2))
#expression(Super^"script text"), round(valorR2,2))
OP <- par(pty="s", mar=c(4,2,4,1))
plot(IDW.out ~ sitiosCOS$p0a5,  main = titleGraph, col.main="blue", asp=3, 
     xlab="Observed values", ylab="Predicted values", pch=16, col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ sitiosCOS$p0a5), col="red", lw=1,lty=1)
abline(0,1)
par(OP)
par


write.csv(rNew@data@values,"Values5-15.csv", row.names = FALSE)

