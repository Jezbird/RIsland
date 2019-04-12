install.packages(c("devtools","Distance","dsm","mads","tidyverse","forcats","knitr","kableExtra"))

install.packages("mrds", type = "source")
library(devtools)
install_github("DistanceDevelopment/Distance")
library(Distance)
library(dsm)
library(mads)
library(tidyverse)
library(forcats)
library(knitr)
library(kableExtra)

?Distance
?mrds
?dsm

## Set working directory
getwd()
setwd(dir = "E:/Users/Jez/Documents/PhD/Data/Island-wide survey/R")
getwd()

prions <- read.csv("AP_brws_multicov.csv")
head(prions)

## First, using the pooled data test different models and adjustments to select the best fit (lowest AIC)

##key key function to use; "hn" gives half-normal (default), "hr" gives hazard-rate and
##"unif" gives uniform. Note that if uniform key is used, covariates cannot be
##included in the model.
##adjustment terms to use; "cos" gives cosine (default), "herm" gives Hermite
##polynomial and "poly" gives simple polynomial. "cos" is recommended. A
##value of NULL indicates that no adjustments are to be ???tted.

## Models to test:
## Half-normal
## Hazard-rate
## Half-normal + Cosine
## Hazard-rate + Cosine
## Uniform + Cosine
## Half-normal + Hermite
## Hazard-rate + Hermite
## Uniform + Hermite
## Half-normal + Polynomial
## Hazard-rate + Polynomial
## Uniform + Polynomial
## Half-normal + Stratum
## Hazard-rate + Stratum
## Half-normal + Observer
## Hazard-rate + Observer
## Half-normal + Veg height
## Hazard-rate + Veg height
## Half-normal + Stratum + Observer
## Hazard-rate + Stratum + Observer
## Half-normal + Stratum + Veg height
## Hazard-rate + Stratum + Veg height
## Half-normal + Observer + Veg height
## Hazard-rate + Observer + Veg height
## Half-normal + Stratum + Observer + Veg height
## Hazard-rate + Stratum + Observer + Veg height

prion.trunc <- 10
prion.half.norm <- ds(prions, truncation=prion.trunc,
                           key="hn",  adjustment=NULL)
prion.haz.rate <- ds(prions, truncation=prion.trunc,
                      key="hr",  adjustment=NULL)
prion.half.norm.cos <- ds(prions, truncation=prion.trunc,
                      key="hn",  adjustment="cos")
prion.haz.rate.cos <- ds(prions, truncation=prion.trunc,
                     key="hr",  adjustment="cos")
prion.unif.cos <- ds(prions, truncation=prion.trunc,
                 key="unif",  adjustment="cos")
prion.half.norm.herm <- ds(prions, truncation=prion.trunc,
                      key="hn",  adjustment="herm")
prion.haz.rate.herm <- ds(prions, truncation=prion.trunc,
                     key="hr",  adjustment="herm")
prion.unif.herm <- ds(prions, truncation=prion.trunc,
                 key="unif",  adjustment="herm")
prion.half.norm.poly <- ds(prions, truncation=prion.trunc,
                      key="hn",  adjustment="poly")
prion.haz.rate.poly <- ds(prions, truncation=prion.trunc,
                     key="hr",  adjustment="poly")
prion.unif.poly <- ds(prions, truncation=prion.trunc,
                 key="unif",  adjustment="poly")

prions$stratum <- ifelse(prions$Region.Label=="1", "1", 
                         (ifelse(prions$Region.Label=="2", "2", 
                                 (ifelse(prions$Region.Label=="3", "3", 
                                         (ifelse(prions$Region.Label=="4", "4",
                                                 (ifelse(prions$Region.Label=="5", "5",
                                                         (ifelse(prions$Region.Label=="6", "6",
                                                                 (ifelse(prions$Region.Label=="7", "7", "8")))))))))))))
prion.strat.hn <- ds(prions, truncation=prion.trunc, quiet=TRUE,
                            formula=~stratum,
                            key="hn")
prion.strat.hr <- ds(prions, truncation=prion.trunc, quiet=TRUE,
                            formula=~stratum,
                            key="hr")
prion.obs.hn <- ds(prions, truncation=prion.trunc,
                      key="hn",  formula=~surveyor)
prion.obs.hr <- ds(prions, truncation=prion.trunc,
                   key="hr",  formula=~surveyor)
prion.veg_height.hn <- ds(prions, truncation=prion.trunc,
                   key="hn",  formula=~veg_height)
prion.veg_height.hr <- ds(prions, truncation=prion.trunc,
                          key="hr",  formula=~veg_height)
prion.strat.obs.hn <- ds(prions, truncation=prion.trunc,
                   key="hn",  formula=~stratum + surveyor)
prion.strat.obs.hr <- ds(prions, truncation=prion.trunc,
                         key="hr",  formula=~stratum + surveyor)
prion.strat.veg.hn <- ds(prions, truncation=prion.trunc,
                         key="hn",  formula=~stratum + veg_height)
prion.strat.veg.hr <- ds(prions, truncation=prion.trunc,
                         key="hr",  formula=~stratum + veg_height)
prion.strat.veg.obs.hn <- ds(prions, truncation=prion.trunc,
                         key="hn",  formula=~stratum + veg_height + surveyor)
prion.strat.veg.obs.hr <- ds(prions, truncation=prion.trunc,
                         key="hr",  formula=~stratum + veg_height + surveyor)

summarize_ds_models(prion.half.norm,prion.haz.rate,prion.half.norm.cos,prion.haz.rate.cos,
                    prion.half.norm.herm,prion.haz.rate.herm,prion.unif.herm,
                    prion.half.norm.poly,prion.haz.rate.poly,prion.unif.poly,prion.strat.hn,
                    prion.strat.hr,prion.obs.hn,prion.obs.hr,prion.veg_height.hn,prion.veg_height.hr,
                    prion.strat.obs.hn,prion.strat.obs.hr,prion.strat.veg.hn,prion.strat.veg.hr,
                    prion.strat.veg.obs.hn,prion.strat.veg.obs.hr, output = "latex")

?summarize_ds_models
?integrate
?distpdf
# Plot simple model
plot(prion.strat.veg.hr, nc = 20, main = "Hazard-rate model with stratum and veg height as covariates",
     pch = 20)

## Now compare all AIC values for the models above. Is there a way of getting R to tabulate all of these?
## The number of adjustments seems to be important - I don't understand what these are - need to find out.

## The Rdistance automated analysis selects a half-normal model with 3 hermite expansions as the best model with
## AIC value of 8779.4773. Don't understand why AIC values are so different when the data is the same. Also why
## the 'best' model is a different one when I use a different R package.


## Now test the effect different covariates have on the data:
## Observer. NB this section of code isn't working. It keeps returning the error: Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
## contrasts can be applied only to factors with 2 or more levels


## Based on the above use a hazard rate model with cosine adjustment to fully stratify the data, 
## calculating separate detection functions for each stratum. 
## Need to figure out justification for truncation - 10 m looks sensible from the graph (see the Rdistance analysis).

prion.full.strat1 <- ds(prions[prions$Region.Label=="1",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")
prion.full.strat2 <- ds(prions[prions$Region.Label=="2",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")
prion.full.strat3 <- ds(prions[prions$Region.Label=="3",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")
prion.full.strat4 <- ds(prions[prions$Region.Label=="4",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")
prion.full.strat5 <- ds(prions[prions$Region.Label=="5",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")
prion.full.strat6 <- ds(prions[prions$Region.Label=="6",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")
prion.full.strat7 <- ds(prions[prions$Region.Label=="7",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")
prion.full.strat8 <- ds(prions[prions$Region.Label=="8",], truncation=prion.trunc, 
                        key="hr",  adjustment="cos")


## The above analysis can't fit detection models for strata 4 or 7 which only had 32
## and 3 detections respectively. Does this affect the AIC value calculation?


## The second of these next two tables generates density estimates for each stratum. Multiply that by the
## area in the first table to get an island total. Results look dodgy. Area 7 which was barely surveyed as
## it is unsuitable cliff habitat had 3 detections resulting in a massively skewed contribution to the
## overall population estimate. The other strata appear to have very low population estimates. There must
## be something weird going on - effort comes out at 261 km but only 157 km of transect were completed. 
## Nowhere has the analysis required the lengths of all the transects where no burrows were detected.

kable(prion.strat.veg.hr$dht$individuals$summary, format = "markdown")

kable(prion.strat.veg.hr$dht$individuals$D,format="markdown")

## Now pool all the data from all strata. 

prion.pooledf0 <- ds(prions, truncation=prion.trunc,
                     key="hr",  adjustment="cos")

kable(prion.pooledf0$dht$individuals$summary, format = "markdown")

kable(prion.pooledf0$dht$individuals$D,format="markdown")

## Compare AIC values from the three approaches
## Separate detection functions overall AIC = 5472.807 - NB this should have AIC values added from the two missing strata
## Stratum as a covariate AIC = 5543.304
## Pool data AIC = 5568.943
## Use the model with stratum as a covariate

## Plot stratum-specific detection functions when using stratum as a covariate. Need to understand how to check 
## goodness of fit of the model.

## Cramer von-Mises test p-value seems to suggest a poor fit and that the data are significantly different from the expected
## Might be worth exploring ditching some strata.

par(mfrow=c(1,2))
plot(prion.strat.veg.hr, main="Prion burrows, \ndetection function uses stratum and veg height as covariates")
covar.fit <- ddf.gof(prion.strat.veg.hr$ddf)
message <- paste("Cramer von-Mises W=", round(covar.fit$dsgof$CvM$W,3), 
                 "\nP=", round(covar.fit$dsgof$CvM$p,3))
text(0.6, 0.1, message, cex=0.8)
