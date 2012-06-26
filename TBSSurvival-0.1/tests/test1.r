# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2012 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## This code is used for testing purposes. The TBSSurvival library does not
## depend on it for any of its functionalities
source("loadlibs.r")
loadlibs()
set.seed(1)

####################
## simple test with the GBSG2 (German Breast Cancer Group 2) data set from the ipred package
library(ipred)
data(GBSG2)
cat('Running MLE on GBSG2 (from ipred package) without covariates\n')
s=tbs.survreg.mle(Surv(GBSG2$time,GBSG2$cens==1) ~ 1,dist="norm",verbose=TRUE) ## try all methods

####################
## test with the colon data set from the survival package
library(survival)
data(colon)
cat('Running MLE on colon (from survival package) without covariates\n')
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ 1,dist="norm",method="Rsolnp",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ 1,dist="norm",method="BFGS",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ 1,dist="norm",method="Nelder-Mead",verbose=TRUE)

## with covariate
cat('Running MLE on colon (from survival package) with covariate=age60\n')
colon$age60=as.numeric(colon$age>60) #threshold defined from medical papers
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",method="Rsolnp",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",method="BFGS",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",method="Nelder-Mead",verbose=TRUE)

s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$node4,dist="norm",method="BFGS",verbose=TRUE)

