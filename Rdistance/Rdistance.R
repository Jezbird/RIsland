
## Set working directory
getwd()
setwd(dir = "C:/Users/Jez/Documents/PhD/Data/Island-wide survey/R/Rdistance")
getwd()

require(Rdistance)

ap.detections <- read.csv("AP_brws.csv")
transects <- read.csv("transect_lengths.csv")

transects[["ufcm"]] = as.factor(transects[["ufcm"]])

## Draw histogram of data - use to set truncate distance
hist(ap.detections$dist, col="grey", main="", xlab="distance (m)")
rug(ap.detections$dist)

## Generate summary stats
summary(ap.detections$dist)

## Automatically detect best distance function
auto <- autoDistSamp(formula = dist~1,
                     detectionData = ap.detections,
                     siteData = transects,
                     area = 10000,
                     likelihoods = c("halfnorm", "hazrate", "uniform", "negexp", "Gamma"),
                     expansions = 0:3,
                     series = c("cosine", "hermite", "simple"),
                     ci = 0.95,
                     R = 500,
                     plot = FALSE,
                     plot.bs = FALSE,
                     w.hi = 150)

## Compare covariates:

## Observer:
hist(ap.detections$dist[transects$observer ==
                                 "JB"], col="grey", main="JB",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                                "JB"],quiet = TRUE)

hist(ap.detections$dist[transects$observer ==
                          "PP"], col="grey", main="PP",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "PP"],quiet = TRUE)

## Veg height
hist(ap.detections$dist[transects$vegheight ==
                          "Short"], col="grey", main="Short",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "Short"],quiet = TRUE)

hist(ap.detections$dist[transects$vegheight ==
                          "Tall"], col="grey", main="Tall",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "Tall"],quiet = TRUE)

## Stratum
hist(ap.detections$dist[transects$ufcm ==
                          "1"], col="grey", main="1",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "1"],quiet = TRUE)

hist(ap.detections$dist[transects$ufcm ==
                          "2"], col="grey", main="2",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "2"],quiet = TRUE)

hist(ap.detections$dist[transects$ufcm ==
                          "3"], col="grey", main="3",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "3"],quiet = TRUE)

hist(ap.detections$dist[transects$ufcm ==
                          "4"], col="grey", main="4",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "4"],quiet = TRUE)

hist(ap.detections$dist[transects$ufcm ==
                          "5"], col="grey", main="5",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "5"],quiet = TRUE)

hist(ap.detections$dist[transects$ufcm ==
                          "6"], col="grey", main="6",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "6"],quiet = TRUE)

hist(ap.detections$dist[transects$ufcm ==
                          "7"], col="grey", main="7",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "7"],quiet = TRUE)

hist(ap.detections$dist[transects$ufcm ==
                          "7"], col="grey", main="7",
     xlab="distance (m)", breaks = 15, xlim = c(0,10))
rug(ap.detections$dist[transects$observer ==
                         "7"],quiet = TRUE)



?hist

autoCov <- autoDistSamp(formula = dist~observer,
                        detectionData = ap.detections,
                        siteData = transects,
                        likelihoods = c("halfnorm", "hazrate", "negexp"),
                        expansions = 0,
                        area = 10000,
                        ci = 0.95,
                        R = 500,
                        plot = FALSE,
                        plot.bs = FALSE,
                        w.hi = 150)


## The code above trials different detection functions and picks the one that best fits your data (lowest AIC value).
## The abundance estimate returned corresponds to a unit area. The area argument is used to multiply this to a new area
## so in this case the original unit of 1m2 is multiplied by 10,000 to give abundance per hectare.
## R is the number of bootstrap iterations used to generate confidence intervals.
## The effective strip width (ESW) is the key information from the detection function that 
## will be used to next estimate abundance (or density). The ESW is calculated by integrating 
## under the detection function. A survey with imperfect detection and ESW equal to X effectively 
## covers the same area as a study with perfect detection out to a distance of X.


## The code below calculates an overall density estimate in burrows per m2  when stratum is included as a covariate
## but it does not calculate stratum-specific density estimates. Rather than by extrapolating overall for the whole
## island area we should extrapolate based upon effort per stratum - need to think about this.

dfuncStratum<- dfuncEstim(formula = dist~ufcm,
                             detectionData = ap.detections,
                             siteData = transects,
                             likelihood = "halfnorm", w.hi = 10)
plot(dfuncStratum, newdata = data.frame(ufcm=
                                             levels(transects$ufcm)), nbins = 10)

covarFit <- abundEstim(dfuncStratum,
                       detectionData=ap.detections,
                       siteData=transects,
                       area=1,
                       ci=0.95,
                       R=20,
                       showProgress = TRUE)
covarFit

