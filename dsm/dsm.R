library(dsm)
library(Distance)
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.8-7. For overview type 'help("mgcv-package")'.
## Loading required package: mrds
## This is mrds 2.1.15
## Built: R 3.2.2; ; 2015-11-24 17:46:38 UTC; unix
## This is dsm 2.2.11
## Built: R 3.2.2; ; 2015-11-24 18:04:45 UTC; unix
library(ggplot2)

# plotting options
gg.opts <- theme(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())

# make the results reproducible
set.seed(11123)

install.packages(c("dsm", "Distance", "knitr", "captioner", "ggplot2", "rgdal",
                     "maptools", "plyr", "tweedie"))

## Set working directory
setwd(dir = "C:/Users/Jez/Documents/PhD/Data/Island-wide survey/R/dsm")
getwd()

segdata <- read.csv("segdata.csv")
head(segdata)

distdata <- read.csv("distdata.csv")

obsdata <- read.csv("obsdata.csv")


## generate detection function

detfc.hr.null<-ds(distdata, max(distdata$distance), key="hr", adjustment=NULL)

## Fitting simple dsm
dsm.xy <- dsm(count~s(x,y), detfc.hr.null, segdata, obsdata, method="REML")

summary(dsm.xy)

vis.gam(dsm.xy, plot.type="contour", view=c("x","y"), asp=1, type="response", contour.col="black", n.grid=100)

dsm.xy.elevation <- dsm(count~s(x,y,k=10) + s(elevation_mean,k=20), detfc.hr.null, segdata, obsdata, group = TRUE, segment.area = "area", method="REML")
summary(dsm.xy.elevation)

plot(dsm.xy.elevation, select=2)

?dsm

dsm.xy.ele.slope <- dsm(count~s(x,y,k=10) + s(elevation_mean,k=20) +s(slope_mean, k=20), detfc.hr.null, segdata, obsdata, group = TRUE, segment.area = "area", method="REML")
summary(dsm.xy.ele.slope)

plot(dsm.xy.ele.slope, select=c(2,3))

gam.check(dsm.xy)
gam.check(dsm.xy.ele.slope)
