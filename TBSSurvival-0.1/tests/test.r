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

w <- options("warn")
options("warn" = -1)
if(require("TBSSurvival",quietly=TRUE)==FALSE) {
  require("survival")
  require("R.utils")
  require("mcmc")
  require("normalp")
  require("eha")
  require("e1071")
  require("coda")
  source("../R/tbs.survreg.bayes.r")
  source("../R/tbs.r")
  source("../R/tbs.survreg.mle.r")
  source("../R/local.r")
}
options("warn" = w[[1]])

library(ipred)
data(GBSG2)
s=tbs.survreg.mle(Surv(GBSG2$time,GBSG2$cens==1) ~ 1,dist="norm",method="Rsolnp",verbose=TRUE)
s=tbs.survreg.mle(Surv(GBSG2$time,GBSG2$cens==1) ~ 1,dist="norm",method="BFGS",verbose=TRUE)
s=tbs.survreg.mle(Surv(GBSG2$time,GBSG2$cens==1) ~ 1,dist="norm",method="Nelder-Mead",verbose=TRUE)

library(survival)
data(colon)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ 1,dist="norm",method="Rsolnp",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ 1,dist="norm",method="BFGS",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ 1,dist="norm",method="Nelder-Mead",verbose=TRUE)
colon$age60=as.numeric(colon$age>60) #from medical paper
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",method="Rsolnp",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",method="BFGS",verbose=TRUE)
s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",method="Nelder-Mead",verbose=TRUE)


######################
# simulation: function used to generate data from Weibull distribution and
#             then compare the Weibull model with TBSmodel.
#             Also, this function is used to analyse if TBSmodel MLE is 
#             working well.
simulation <- function(gen=list(n=100,cens=0,shape=2,scale=3),
                       est=list(dist="norm"),n.copies=1,initial.seed=1234,method="Rsonlp") {
  # gen: is a list with the parameters to generate from a Weibull distribution
  #          n: sample size
  #       cens: censor rate
  #      shape: shape Weibull parameter
  #      scale: scale Weibull parameter
  # est: is a list with the parameters used to estimate the model
  #      dist: error distribution of TBS model
  # n.copies: how many copies will be performed in the simulation
  # initial.seed: to control the simulation

  result.1  <- matrix(est$dist,n.copies,1)
  result.2  <- matrix(NA,n.copies,20)

  for (i in 1:n.copies) {
    # Generating data
    initial.time <- .gettime()
    seed <- initial.seed+i-1
    set.seed(seed)
    d <- list()
    d$time  <- rweibull(gen$n,shape=1/gen$shape,scale=exp(gen$scale))
    d$delta <- rep(1,gen$n)
    censor  <- quantile(d$time,probs=(1-gen$cens))
    d$delta <- ifelse(d$time > censor, 0, 1)
    d$time  <- ifelse(d$time > censor, censor, d$time)

    # Original Survival function
    orig.y   <- pweibull(d$time[d$delta == 1],shape=1/gen$shape,scale=exp(gen$scale))

    # Weibull fit
    weib.fit <- survreg(Surv(d$time,d$delta) ~ 1, dist="weibull")
    weib.y   <- pweibull(d$time[d$delta == 1],shape=1/weib.fit$scale,
                                              scale=exp(weib.fit$coefficients))

    # Tbs fit
    tbs.fit  <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1,dist=est$dist,method=method)

    # Result matrix
    if (tbs.fit$convergence) {
      tbs.y <- ptbs(d$time[d$delta == 1],lambda=tbs.fit$par[1],xi=tbs.fit$par[2],
                    beta=tbs.fit$par[3],dist=est$dist)
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$shape,gen$scale,tbs.fit$par,
                        max(abs(orig.y-tbs.y)),mean(abs(orig.y-tbs.y)),
                        1/weib.fit$scale,exp(weib.fit$coefficients),
                        max(abs(orig.y-weib.y)),mean(abs(orig.y-weib.y)),
                        ifelse(max(abs(orig.y-tbs.y)) > max(abs(orig.y-weib.y)), 0, 1),
                        ifelse(mean(abs(orig.y-tbs.y)) > mean(abs(orig.y-weib.y)), 0, 1),
                        gen$scale*log(2)^(1/gen$shape),
                        exp(tbs.fit$par[3]),
                        exp(weib.fit$coefficients)*log(2)^weib.fit$scale,run.time)
    } else {
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$shape,gen$scale,rep(NA,5),
                        1/weib.fit$scale,exp(weib.fit$coefficients),
                        max(abs(orig.y-weib.y)),mean(abs(orig.y-weib.y)),
                        rep(NA,2),gen$scale*log(2)^(1/gen$shape),NA,
                        exp(weib.fit$coefficients)*log(2)^weib.fit$scale,run.time)
    }
  }
  
  colnames(result.1) <- "dist"
  colnames(result.2) <- c("gen.seed","n","censor","orig.par1","orig.par2",
                          "tbs.lambda","tbs.xi","tbs.beta","tbs.Max.AE",
                          "tbs.Mean.AE","weib.par1","weib.par2",
                          "weib.Max.AE","weib.Mean.AE","comp.Max.AE","comp.Mean.AE",
                          "orig.median","tbs.median","weib.median","run.time")

  return(result.2)
}


simulation(gen=list(n=1000,cens=0,scale=2,shape=3),est=list(dist="norm"),
           n.copies=1,initial.seed=1235,method="BFGS")


