#Proyecto final: Estimación del riesgo relativo de SIDS utilizando un modelo CAR

library(tidyverse); library(rgdal); library(maptools); library(maps); library(spdep); 
library(tmap); library(GISTools); library(sp); library(RColorBrewer); library(sf);
library(graphics); library(raster); library(DCluster); library(spatialreg); library(leaflet);
library(ggpubr); library(epitools); library(gplots); library(mgcv)
library(dplyr); library(ggplot2); library(scales); library(magrittr); library(gstat);
library(gridExtra); library(ggthemes); library(tidyverse); library(maptools); library(rgdal); 
library(lattice); library(RSAGA); library(spatstat); library(cowplot); library(geoR); library(kableExtra);
library(stargazer); library(gamlss.spatial)

#Se importa y maneja la base
setwd("C:/Users/verdu/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESPECIALIDAD/SEGUNDO SEMESTRE/ESPACIAL/sids")

nc <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
st_crs(nc) <- "+proj=longlat +datum=NAD27"
row.names(nc) <- as.character(nc$FIPSNO)
sf_use_s2(FALSE)
plot(st_geometry(nc), axes=TRUE)
text(st_coordinates(st_centroid(st_geometry(nc), of_largest_polygon=TRUE)), label=nc$FIPSNO, cex=0.5)

t1<-nc%>%as.data.frame%>%dplyr::select(BIR74, SID74, BIR79, SID79)%>%finalfit::ff_glimpse()
kable(t1, format="latex")

gal_file <- system.file("weights/ncCR85.gal", package="spData")[1]
ncCR85 <- read.gal(gal_file, region.id=nc$FIPSNO)
ncCR85

gal_file <- system.file("weights/ncCC89.gal", package="spData")[1]
ncCC89 <- read.gal(gal_file, region.id=nc$FIPSNO)
ncCC89

r.id <- attr(ncCC89, "region.id")
ncCC89[[match("37001", r.id)]]
r.id[ncCC89[[match("37001", r.id)]]]
as.character(nc$NAME)[card(ncCC89) == 0]

#### Estimaciones para 1974
#Cálculo de SIRS para 1974 (SMR por standardized mortality ratio)
nc$Observed <- nc$SID74
nc$Population <- nc$BIR74
r <- sum(nc$Observed)/sum(nc$Population)
nc$Expected <- nc$Population * r
nc$SMR <- nc$Observed/nc$Expected

#Gráfica de intervalos de confianza exactos para los SMR por método Poisson
ggplot(nc, aes(y = SMR, x = NAME)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymax = pois.exact(nc$SMR)$upper, ymin = pois.exact(nc$SMR)$lower))+
  geom_hline(yintercept=1, aes(color="red"))+
  theme_minimal()+
  xlab("Condado")+
  ylab("Riesgo relativo")+
  theme(axis.text.x = element_text(angle = 90))


#Observamos las características de las vecindades y el índice de Moran
col.W <- nb2listw(ncCR85, zero.policy = TRUE)
moranI.test(Observed ~ offset(log(Expected)), as.data.frame(nc), "negbin", 999, listw = col.W, n = length(ncCR85), S0 = Szero(col.W))
plot(st_geometry(nc), border="grey")
plot(ncCR85, st_centroid(st_geometry(nc), of_largest_polygon), add=TRUE, col="blue")
#Se confirma la autocorrelación espacial a través de la prueba de Moran

#Calculamos los modelos
##1. Modelo Poisson-Gamma
eb <- empbaysmooth(nc$Observed, nc$Expected)
nc$EBPG <- eb$smthrr

##2. Modelo Log-normal
ebln <- lognormalEB(nc$Observed, nc$Expected)
nc$EBLN <- exp(ebln$smthrr)

## Modelo SAR
sar.nc<-lagsarlm(SMR~1, data=nc, listw=col.W)
summary(sar.nc)

nc$SAR<-sar.nc$fitted.values
nc$SAR.res<-sar.nc$residuals

## Modelo por Campo Aleatorio de Markov
nc_poly1<-nb2nb(poly2nb(nc))
nc$r.id<-as.factor(r.id)

## Primero un MRF de rango completo simple
b <- gam(SMR ~ s(r.id, bs="mrf",k=20,xt=nc_poly1), data=nc, method="REML", family=nb())
nc$mrf<-b$fitted.values

## Gráficos
qpal<-colorNumeric("RdBu", nc$Observed) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal(Observed)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal, values = ~nc$Observed,
            title = "Valores Observados",
            opacity = 1
  )

qpal2<-colorNumeric("RdBu", nc$Expected) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal2(Expected)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal2, values = ~nc$Expected,
            title = "Valores Esperados",
            opacity = 1
  )

qpal3<-colorNumeric("RdBu", nc$EBLN) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal3(EBLN)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal3, values = ~nc$EBLN,
            title = "Poisson-Gamma",
            opacity = 1
  )

qpal4<-colorNumeric("RdBu", nc$EBPG) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal4(EBPG)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal4, values = ~nc$EBPG,
            title = "LogNormal",
            opacity = 1
  )

qpal5<-colorNumeric("RdBu", nc$SAR) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal5(SAR)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal5, values = ~nc$SAR,
            title = "SAR",
            opacity = 1
  )

qpal6<-colorNumeric("RdBu", nc$SAR.res) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal6(SAR.res)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal6, values = ~nc$SAR.res,
            title = "SAR residuals",
            opacity = 1
  )

qpalb<-colorNumeric("RdBu", nc$mrf) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpalb(mrf)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpalb, values = ~nc$mrf,
            title = "MRF",
            opacity = 1
  )

#### Estimaciones para 1978

#Cálculo de SIRS para 1978 (SMR por standardized mortality ratio)
nc$Observed79 <- nc$SID79
nc$Population79 <- nc$BIR79
r79 <- sum(nc$Observed79)/sum(nc$Population79)
nc$Expected79 <- nc$Population79 * r79
nc$SMR79 <- nc$Observed79/nc$Expected79

#Gráfica de intervalos de confianza exactos para los SMR por método Poisson
ggplot(nc, aes(y = SMR79, x = NAME)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymax = pois.exact(nc$SMR79)$upper, ymin = pois.exact(nc$SMR79)$lower))+
  geom_hline(yintercept=1, aes(color="red"))+
  theme_minimal()+
  xlab("Condado")+
  ylab("Riesgo relativo")+
  theme(axis.text.x = element_text(angle = 90))


#Observamos las características de las vecindades y el índice de Moran

moranI.test(Observed79 ~ offset(log(Expected79)), as.data.frame(nc), "negbin", 999, listw = col.W, n = length(ncCR85), S0 = Szero(col.W))
#Se confirma la autocorrelación espacial a través de la prueba de Moran

#Calculamos los modelos
##1. Modelo Poisson-Gamma
eb79 <- empbaysmooth(nc$Observed79, nc$Expected79)
nc$EBPG79 <- eb79$smthrr

##2. Modelo Log-normal
ebln79 <- lognormalEB(nc$Observed79, nc$Expected79)
nc$EBLN79 <- exp(ebln79$smthrr)

## Modelo SAR
sar.nc79<-lagsarlm(SMR79~1, data=nc, listw=col.W)
summary(sar.nc79)

nc$SAR79<-sar.nc79$fitted.values
nc$SAR.res79<-sar.nc79$residuals

## Modelo por Campo Aleatorio de Markov
## Primero un MRF de rango completo simple
b79 <- gam(SMR79 ~ s(r.id, bs="mrf",xt=nc_poly1), data=nc, method="REML", family=nb())
nc$mrf79<-b79$fitted.values

plot(b79$fitted.values,scheme=1)

##Gráficos

qpal6<-colorNumeric("PRGn", nc$Observed79) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal6(Observed79)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal6, values = ~nc$Observed79,
            title = "Valores Observados 1979",
            opacity = 1
  )

qpal7<-colorNumeric("PRGn", nc$Expected79) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal7(Expected79)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal7, values = ~nc$Expected79,
            title = "Valores Esperados 1979",
            opacity = 1
  )

qpal8<-colorNumeric("PRGn", nc$EBLN79) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal8(EBLN79)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal8, values = ~nc$EBLN79,
            title = "Poisson-Gamma 1979",
            opacity = 1
  )

qpal9<-colorNumeric("PRGn", nc$EBPG79) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal9(EBPG79)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal9, values = ~nc$EBPG79,
            title = "LogNormal",
            opacity = 1
  )

qpal10<-colorNumeric("PRGn", nc$SAR79) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal10(SAR79)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal10, values = ~nc$SAR79,
            title = "SAR",
            opacity = 1
  )

qpal11<-colorNumeric("PRGn", nc$SAR.res79) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpal11(SAR.res79)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpal11, values = ~nc$SAR.res79,
            title = "SAR residuals",
            opacity = 1
  )

qpalb79<-colorNumeric("PRGn", nc$mrf79) 
leaflet(nc) %>%
  addPolygons(stroke = FALSE, fillOpacity = .8, smoothFactor = 0.2, color = ~qpalb79(mrf79)) %>%
  addTiles() %>%
  addLegend("bottomright", pal = qpalb79, values = ~nc$mrf79,
            title = "MRF",
            opacity = 1
  )

AIC(b)
AIC(sar.nc)

AIC(b79)
AIC(sar.nc79)

ggplot(data=nc, aes(x=BIR74, color="blue"))+
  geom_histogram(aes(x=SID74, color="red", position="dodge"))+
  theme_minimal()

plot(nc[,c(9, 10, 12, 13, 27, 37)])
legend("bottomright", legend = "Escala", pch = 21)

