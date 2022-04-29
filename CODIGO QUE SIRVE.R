library(dplyr); library(ggplot2); library(scales); library(magrittr); library(sp); library(gstat); library(gridExtra)
library(ggthemes); library(tidyverse); library(maptools); library(rgdal); library(lattice); library(RSAGA); library(spatstat)
library(ggpubr); library(cowplot); library(raster); library(sf); library(geoR); library(kableExtra);
library(stargazer)

setwd("C:/Users/verdu/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESPECIALIDAD/SEGUNDO SEMESTRE/ESPACIAL/INTENTO 2")

#Importando el dataframe

data(meuse)

#ACA INSERTAR DESCRIPTIVOS
t1<-meuse%>%as.data.frame%>%finalfit::ff_glimpse()
kable(t1, format="latex")

#'attach coordinates' - convert table to a point map:
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")

#Visualizamos río
data(meuse.riv)
tmp <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
meuse.riv <- SpatialPolygons(tmp)
zinc.plt<-spplot(meuse, "zinc", scales=list(draw=T), do.log=TRUE, colorkey=TRUE, sp.layout=list("sp.polygons", meuse.riv, col="blue"),
       main="Concentración de Zinc (ppm)")
cadmium.plt<-spplot(meuse, "cadmium", scales=list(draw=T), do.log=TRUE, colorkey=TRUE, sp.layout=list("sp.polygons", meuse.riv, col="blue"),
       main="Concentración de Cadmio (ppm)")
lead.plt<-spplot(meuse, "lead", scales=list(draw=T), do.log=TRUE, colorkey=TRUE, sp.layout=list("sp.polygons", meuse.riv, col="blue"),
       main="Concentración de Plomo (ppm)")
copper.plt<-spplot(meuse, "copper", scales=list(draw=T), do.log=TRUE, colorkey=TRUE, sp.layout=list("sp.polygons", meuse.riv, col="blue"),
       main="Concentración de Cobre (ppm)")
print(zinc.plt, split=c(1, 1, 2, 2), more=TRUE)
print(cadmium.plt, split=c(2, 1, 2, 2), more=TRUE)
print(lead.plt, split=c(1, 2, 2, 2), more=TRUE)
print(copper.plt, split=c(2, 2, 2, 2), more=TRUE)



download.file("http://spatial-analyst.net/book/system/files/meuse.zip", destfile=paste(getwd(), "meuse.zip", sep="/"))
grid.list <- c("ahn.asc", "dist.asc", "ffreq.asc", "soil.asc")

# unzip the maps in a loop:
for(j in grid.list){
  fname <- unzip(files=j, zipfile="meuse.zip")
}


meuse.grid <- readGDAL(grid.list[1])
names(meuse.grid)[1] <- sub(".asc", "", grid.list[1])
for(i in grid.list[-1]) {
  meuse.grid@data[sub(".asc", "", i[1])] <- readGDAL(paste(i))$band1
}

#Ponemos el sistema de coordenadas correcto
#Contamos como factores a las variables categóricas
meuse.grid$ffreq <- as.factor(meuse.grid$ffreq)
meuse.grid$soil <- as.factor(meuse.grid$soil)

ffreq.plt <- spplot(meuse.grid["ffreq"],  col.regions=grey(runif(length(levels(meuse.grid$ffreq)))), main="Frecuencia de inundación")
dist.plt <- spplot(meuse.grid["dist"], col.regions=grey(rev(seq(0,1,0.025))), main="Distancia al río")
soil.plt <- spplot(meuse.grid["soil"],  col.regions=grey(runif(length(levels(meuse.grid$soil)))), main="Tipos de suelo")
ahn.plt <- spplot(meuse.grid["ahn"],  col.regions=grey(rev(seq(0,1,0.025))), main="Elevación [cm]")

print(ffreq.plt, split=c(1, 1, 2, 2), more=TRUE)
print(dist.plt, split=c(2, 1, 2, 2), more=TRUE)
print(ahn.plt, split=c(1, 2, 2, 2), more=TRUE)
print(soil.plt, split=c(2, 2, 2, 2), more=TRUE)

spplot(meuse, "zinc", do.log = T, colorkey = TRUE)

proj4string(meuse.grid) <- CRS("+init=epsg:28992")

meuse@data<-meuse@data%>%select(1,2,3,4,6,7,8,9,10,11,12)
meuse.ov <- over(meuse, meuse.grid)
meuse.ov <- cbind(meuse.ov, meuse[c("zinc", "cadmium","copper","lead")]@data, meuse@coords)

#De acá en adelante, hacer x4

par(mfrow = c(4, 2))
scatter.smooth(meuse.ov$dist, meuse.ov$zinc, span=18/19, col="grey", xlab="Distancia al río (log)", ylab="Zinc (ppm)")
scatter.smooth(meuse.ov$ahn, meuse.ov$zinc, span=18/19,col="grey", xlab="Elevación (cm)", ylab="Zinc (ppm)")

scatter.smooth(meuse.ov$dist, meuse.ov$cadmium, span=18/19, col="blue", xlab="Distancia al río (log)", ylab="Cadmio (ppm)")
scatter.smooth(meuse.ov$ahn, meuse.ov$cadmium, span=18/19,col="blue", xlab="Elevación (cm)", ylab="Cadmio (ppm)")

scatter.smooth(meuse.ov$dist, meuse.ov$copper, span=18/19, col="red", xlab="Distancia al río (log)", ylab="Cobre (ppm)")
scatter.smooth(meuse.ov$ahn, meuse.ov$copper, span=18/19,col="red", xlab="Elevación (cm)", ylab="Cobre (ppm)")

scatter.smooth(meuse.ov$dist, meuse.ov$lead, span=18/19, col="purple", xlab="Distancia al río (log)", ylab="Plomo (ppm)")
scatter.smooth(meuse.ov$ahn, meuse.ov$lead, span=18/19,col="purple", xlab="Elevación (cm)", ylab="Plomo (ppm)")
dev.off()

par(mfrow=c(4,2))
boxplot(log1p(meuse.ov$zinc) ~ meuse.ov$soil, col=grey(runif(length(levels(meuse.ov$soil)))), xlab="Tipos de suelo", ylab="Zinc (ppm)")
boxplot(log1p(meuse.ov$zinc) ~ meuse.ov$ffreq, col=grey(runif(length(levels(meuse.ov$ffreq)))),  xlab="Frecuencia de inundación", ylab="Zinc (ppm)")

boxplot(log1p(meuse.ov$cadmium) ~ meuse.ov$soil, col=grey(runif(length(levels(meuse.ov$soil)))), xlab="Tipos de suelo", ylab="Cadmio (ppm)")
boxplot(log1p(meuse.ov$cadmium) ~ meuse.ov$ffreq, col=grey(runif(length(levels(meuse.ov$ffreq)))),  xlab="Frecuencia de inundación", ylab="Cadmio (ppm)")

boxplot(log1p(meuse.ov$copper) ~ meuse.ov$soil, col=grey(runif(length(levels(meuse.ov$soil)))), xlab="Tipos de suelo", ylab="Cobre (ppm)")
boxplot(log1p(meuse.ov$copper) ~ meuse.ov$ffreq, col=grey(runif(length(levels(meuse.ov$ffreq)))),  xlab="Frecuencia de inundación", ylab="Cobre (ppm)")

boxplot(log1p(meuse.ov$lead) ~ meuse.ov$soil, col=grey(runif(length(levels(meuse.ov$soil)))), xlab="Tipos de suelo", ylab="Plomo (ppm)")
boxplot(log1p(meuse.ov$lead) ~ meuse.ov$ffreq, col=grey(runif(length(levels(meuse.ov$ffreq)))),  xlab="Frecuencia de inundación", ylab="Plomo (ppm)")
dev.off()

pairs(zinc ~ ahn+dist+lead+copper+cadmium, meuse.ov, labels=c("Zinc", "Elevación", "Distancia\nal\nrío",
                                                              "Plomo", "Cobre", "Cadmio"))

cor(meuse.grid$ahn, meuse.grid$dist, use="complete.obs")


lm.zinc <- lm(log1p(zinc) ~ ahn+dist+ffreq+soil, meuse.ov)
summary(lm.zinc)
par(mfrow=c(2,2))
plot(lm.zinc)

lm.cadmium <- lm(log1p(cadmium) ~ ahn+dist+ffreq+soil, meuse.ov)
summary(lm.cadmium)
par(mfrow=c(2,2))
plot(lm.cadmium)

lm.copper <- lm(log1p(copper) ~ ahn+dist+ffreq+soil, meuse.ov)
summary(lm.copper)
par(mfrow=c(2,2))
plot(lm.copper)

lm.lead <- lm(log1p(lead) ~ ahn+dist+ffreq+soil, meuse.ov)
summary(lm.lead)
par(mfrow=c(2,2))
plot(lm.lead)

stargazer(lm.zinc, lm.cadmium, lm.copper, lm.lead, ci=T)


coordinates(meuse.ov)<-~x+y
proj4string(meuse.ov) <- CRS("+init=epsg:28992")


####ZINC####
zinc.svar <- variogram(log1p(zinc) ~ 1, meuse)

zinc.rsvar <- variogram(residuals(lm.zinc) ~ 1, data=meuse.ov)
zinc.ivgm <- vgm(nugget=0, model="Exp", range=sqrt(diff(meuse@bbox["x",])^2 + diff(meuse@bbox["y",])^2)/4, psill=var(residuals(lm.zinc)))
zinc.rvgm <- fit.variogram(zinc.rsvar, model=zinc.ivgm)
plot(zinc.rsvar, zinc.rvgm, ylab = "Semivarianza",xlab="Distancia", main="Residuales")

zinc.vgm <- fit.variogram(zinc.svar, model=zinc.ivgm)
plot(zinc.svar, zinc.vgm, pch="+", pl=TRUE, col="black", main="log1p(zinc)", xlab="Distancia", ylab="Semivarianza")
abline(h=0.518, col="red", lty=2, lwd=3)

zinc.rvgm.plt <- plot(zinc.rsvar, zinc.rvgm, pc="+", pl=FALSE, col="black", main="Residuales")
abline(h=0.518, col="red", lty=2, lwd=3)

plot_variogram <- function(v, m) {
  preds = variogramLine(m, maxdist = max(v$dist))
  ggplot() + 
    geom_point(data = v, aes(x = dist, y = gamma, size=np)) +
    geom_line(data = preds, aes(x = dist, y = gamma))+
    theme_hc()
}

zincplot1<-plot_variogram(zinc.rsvar,zinc.rvgm)+
  geom_hline(yintercept=var(log1p(meuse$zinc)), color="red")+
  annotate(geom="text", y=0.55, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

zincplot2<-plot_variogram(zinc.svar,zinc.vgm)+
  geom_hline(yintercept=var(log1p(meuse$zinc)), color="red")+
  annotate(geom="text", y=0.55, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

ggarrange(zincplot2, zincplot1, common.legend = T, labels=c("A: Ordinario", "B: Residuales"))

# synchronize the two plots:
zinc.rvgm.plt$x.limits <- zinc.vgm.plt$x.limits
zinc.rvgm.plt$y.limits <- zinc.vgm.plt$y.limits

print(zinc.vgm.plt, split=c(1,1,2,1), more=TRUE)

print(zinc.rvgm.plt, split=c(2,1,2,1), more=FALSE)

plot(zinc.rvgm.plt)

zinc.rk <- krige(formula=log1p(zinc)~ahn+dist+ffreq+soil, meuse.ov, meuse.grid, zinc.rvgm)
zinc.rk$var1.rk <- expm1(zinc.rk$var1.pred)

at.zinc <- seq(min(meuse$zinc),max(meuse$zinc),sd(meuse$zinc)/5)
zinc.rk.plt <- spplot(zinc.rk["var1.rk"], col.regions=grey(rev(seq(0,0.97,1/length(at.zinc)))), at=at.zinc,
                      main="Zinc", sp.layout=list("sp.points",pch="+", col="black", meuse))
print(zinc.rk.plt)

cross.zinc.rk <- krige.cv(log(zinc)~dist+ffreq+soil+ahn, meuse.ov, zinc.rvgm, verbose=FALSE)
1-var(cross.zinc.rk$residual, na.rm=T)/var(log1p(meuse$zinc))

####CADMIO####
cadmium.svar <- variogram(log1p(cadmium) ~ 1, meuse)

cadmium.rsvar <- variogram(residuals(lm.cadmium) ~ 1, data=meuse.ov)
cadmium.ivgm <- vgm(nugget=0, model="Exp", range=sqrt(diff(meuse@bbox["x",])^2 + diff(meuse@bbox["y",])^2)/4, psill=var(residuals(lm.cadmium)))
cadmium.rvgm <- fit.variogram(cadmium.rsvar, model=cadmium.ivgm)
cadmium.rvgm

cadmium.vgm <- fit.variogram(cadmium.svar, model=cadmium.ivgm)
cadmium.vgm.plt <- plot(cadmium.svar, cadmium.vgm, pch="+", pl=TRUE, col="black", main="log1p(cadmium)")

cadmium.rvgm.plt <- plot(cadmium.rsvar, cadmium.rvgm, pc="+", pl=FALSE, col="black", main="Residuals")

# synchronize the two plots:
cadmium.rvgm.plt$x.limits <- cadmium.vgm.plt$x.limits
cadmium.rvgm.plt$y.limits <- cadmium.vgm.plt$y.limits
print(cadmium.vgm.plt, split=c(1,1,2,1), more=TRUE)
print(cadmium.rvgm.plt, split=c(2,1,2,1), more=FALSE)


cadmium.rk <- krige(formula=log1p(cadmium)~ahn+dist+ffreq+soil, meuse.ov, meuse.grid, cadmium.rvgm)
cadmium.rk$var1.rk <- expm1(cadmium.rk$var1.pred)

at.cadmium <- seq(min(meuse$cadmium),max(meuse$cadmium),sd(meuse$cadmium)/5)
cadmium.rk.plt <- spplot(cadmium.rk["var1.rk"], col.regions=grey(rev(seq(0,0.97,1/length(at.cadmium)))), at=at.cadmium,
                      main="Cadmio", sp.layout=list("sp.points", pch="+", col="black", meuse))
print(cadmium.rk.plt)

cross.cadmium.rk <- krige.cv(log1p(cadmium)~dist+ffreq+soil+ahn, meuse.ov, cadmium.rvgm, verbose=FALSE)
1-var(cross.cadmium.rk$residual, na.rm=T)/var(log1p(meuse$cadmium))

cadplot1<-plot_variogram(cadmium.rsvar,cadmium.rvgm)+
  geom_hline(yintercept=var(log1p(meuse$cadmium)), color="red")+
  annotate(geom="text", y=0.54, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

cadplot2<-plot_variogram(cadmium.svar,cadmium.vgm)+
  geom_hline(yintercept=var(log1p(meuse$cadmium)), color="red")+
  annotate(geom="text", y=0.54, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

ggarrange(cadplot2, cadplot1, common.legend = T, labels=c("A: Ordinario", "B: Residuales"))


####COBRE####
copper.svar <- variogram(log1p(copper) ~1, meuse)

copper.rsvar <- variogram(residuals(lm.copper) ~1, data=meuse.ov)
copper.ivgm <- vgm(nugget=0, model="Exp", range=sqrt(diff(meuse@bbox["x",])^2 + diff(meuse@bbox["y",])^2)/4, psill=var(residuals(lm.copper)))
copper.rvgm <- fit.variogram(copper.rsvar, model=copper.ivgm)
copper.rvgm

copper.vgm <- fit.variogram(copper.svar, model=copper.ivgm)
copper.vgm.plt <- plot(copper.svar, copper.vgm, pch="+", pl=TRUE, col="black", main="log1p(copper)")

copper.rvgm.plt <- plot(copper.rsvar, copper.rvgm, pc="+", pl=FALSE, col="black", main="Residuals")


copper.rvgm.plt$x.limits <- copper.vgm.plt$x.limits
copper.rvgm.plt$y.limits <- copper.vgm.plt$y.limits
print(copper.vgm.plt, split=c(1,1,2,1), more=TRUE)
print(copper.rvgm.plt, split=c(2,1,2,1), more=FALSE)


copper.rk <- krige(formula=log1p(copper)~ahn+dist+ffreq+soil, meuse.ov, meuse.grid, copper.rvgm)
copper.rk$var1.rk <- expm1(copper.rk$var1.pred)

at.copper <- seq(min(meuse$copper),max(meuse$copper),sd(meuse$copper)/5)
copper.rk.plt <- spplot(copper.rk["var1.rk"], col.regions=grey(rev(seq(0,0.97,1/length(at.copper)))), at=at.copper,
                         main="Cobre", sp.layout=list("sp.points", pch="+", col="black", meuse))
print(copper.rk.plt)

cross.copper.rk <- krige.cv(log1p(copper)~dist+ffreq+soil+ahn, meuse.ov, copper.rvgm, verbose=FALSE)
1-var(cross.copper.rk$residual, na.rm=T)/var(log1p(meuse$copper))

copplot1<-plot_variogram(copper.rsvar,copper.rvgm)+
  geom_hline(yintercept=var(log1p(meuse$copper)), color="red")+
  annotate(geom="text", y=0.26, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

copplot2<-plot_variogram(copper.svar,copper.vgm)+
  geom_hline(yintercept=var(log1p(meuse$copper)), color="red")+
  annotate(geom="text", y=0.26, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

ggarrange(copplot2, copplot1, common.legend = T, labels=c("A: Ordinario", "B: Residuales"))

####PLOMO####
lead.svar <- variogram(log1p(lead) ~ 1, meuse)

lead.rsvar <- variogram(residuals(lm.lead) ~ 1, data=meuse.ov)
lead.ivgm <- vgm(nugget=0, model="Exp", range=sqrt(diff(meuse@bbox["x",])^2 + diff(meuse@bbox["y",])^2)/4, psill=var(residuals(lm.lead)))
lead.rvgm <- fit.variogram(lead.rsvar, model=lead.ivgm)

lead.vgm <- fit.variogram(lead.svar, model=lead.ivgm)
lead.vgm.plt <- plot(lead.svar, lead.vgm, pch="+", pl=TRUE, col="black", main="log1p(lead)")

lead.rvgm.plt <- plot(lead.rsvar, lead.rvgm, pc="+", pl=FALSE, col="black", main="Residuals")

lead.rvgm.plt$x.limits <- lead.vgm.plt$x.limits
lead.rvgm.plt$y.limits <- lead.vgm.plt$y.limits
print(lead.vgm.plt, split=c(1,1,2,1), more=TRUE)
print(lead.rvgm.plt, split=c(2,1,2,1), more=FALSE)

lead.rk <- krige(formula=log1p(lead)~ahn+dist+ffreq+soil, meuse.ov, meuse.grid, lead.rvgm)
lead.rk$var1.rk <- expm1(lead.rk$var1.pred)

at.lead <- seq(min(meuse$lead),max(meuse$lead),sd(meuse$lead)/5)
lead.rk.plt <- spplot(lead.rk["var1.rk"], col.regions=grey(rev(seq(0,0.97,1/length(at.lead)))), at=at.lead,
                        main="Plomo", sp.layout=list("sp.points", pch="+", col="black", meuse))
print(lead.rk.plt)

cross.lead.rk <- krige.cv(log1p(lead)~dist+ffreq+soil+ahn, meuse.ov, lead.rvgm, verbose=FALSE)
1-var(cross.lead.rk$residual, na.rm=T)/var(log1p(meuse$lead))

leadplot1<-plot_variogram(lead.rsvar,lead.rvgm)+
  geom_hline(yintercept=var(log1p(meuse$lead)), color="red")+
  annotate(geom="text", y=0.45, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

leadplot2<-plot_variogram(lead.svar,lead.vgm)+
  geom_hline(yintercept=var(log1p(meuse$lead)), color="red")+
  annotate(geom="text", y=0.45, x=200, label="Varianza muestral",
           color="red")+
  xlab("Distancia")+
  ylab("Semivarianza")+
  labs(size="Número de pares")

ggarrange(leadplot2, leadplot1, common.legend = T, labels=c("A: Ordinario", "B: Residuales"))

print(zinc.rk.plt, split=c(1, 1, 2, 2), more=TRUE)
print(cadmium.rk.plt, split=c(2, 1, 2, 2), more=TRUE)
print(lead.rk.plt, split=c(1, 2, 2, 2), more=TRUE)
print(copper.rk.plt, split=c(2, 2, 2, 2), more=TRUE)

variogramas<-ggarrange(zincplot2, zincplot1,
                       cadplot2, cadplot1, 
                       copplot2, copplot1,
                       leadplot2, leadplot1,
                       common.legend = T, 
                       labels=c("A.", "B.",
                                "C.", "D.",
                                "E.", "F.",
                                "G.", "H."), ncol=2, nrow=4)
print(variogramas)

ggsave("variogramas.jpg", plot=variogramas, dpi=300, units="cm", width=30, height=40)
