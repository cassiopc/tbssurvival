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
source('loadlibs.r')
loadlibs()

######################
# simulation: function used to generate data from Weibull distribution and
#             then compare the Weibull model with TBSmodel.
#             Also, this function is used to analyse if TBSmodel MLE is 
#             working well.
sim.weib <- function(gen=list(n=100,cens=0,shape=2,scale=3),est=list(dist="norm"),prefix="",
                     n.copies=1,initial.seed=1234,method="Rsolnp",verbose=TRUE) {
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
  result.2  <- matrix(NA,n.copies,23)

  for (i in 1:n.copies) {
    # Generating data
    initial.time <- .gettime()
    seed <- initial.seed+i-1
    set.seed(seed)
    d <- NULL
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
    tbs.fit  <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1,dist=est$dist,method=method,verbose=verbose)

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
                        exp(gen$scale)*log(2)^gen$shape,median(d$time),
                        exp(tbs.fit$par[3]),
                        (exp(tbs.fit$par[3])-exp(gen$scale)*log(2)^gen$shape)^2,
                        exp(weib.fit$coefficients)*log(2)^weib.fit$scale,
                        (exp(weib.fit$coefficients)*log(2)^weib.fit$scale-
                         exp(gen$scale)*log(2)^gen$shape)^2,run.time)
    } else {
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$shape,gen$scale,rep(NA,5),
                        1/weib.fit$scale,exp(weib.fit$coefficients),
                        max(abs(orig.y-weib.y)),mean(abs(orig.y-weib.y)),
                        rep(NA,2),exp(gen$scale)*log(2)^gen$shape,median(d$time),NA,NA,
                        exp(weib.fit$coefficients)*log(2)^weib.fit$scale,
                        (exp(weib.fit$coefficients)*log(2)^weib.fit$scale-
                         exp(gen$scale)*log(2)^gen$shape)^2,run.time)
    }
    rm(d)
  }
  
  colnames(result.1) <- "dist"
  colnames(result.2) <- c("gen.seed","n","censor","orig.par1","orig.par2",
                          "tbs.lambda","tbs.xi","tbs.beta","tbs.Max.AE",
                          "tbs.Mean.AE","weib.par1","weib.par2",
                          "weib.Max.AE","weib.Mean.AE","comp.Max.AE","comp.Mean.AE",
                          "orig.median","data.median","tbs.median","tbs.sqe",
                          "weib.median","weib.sqe","run.time")

  out <- data.frame(result.1,result.2)
  name <- paste(prefix,"Sim_Weib-",est$dist,"_",gen$n,"_",gen$cens,"_",gen$shape,"_",gen$scale,".csv",sep="")
  write.csv(out,file=name)
#  return(out)
}

######################
# sim.tbs: function used to generate data from TBS model and 
#          then check the estimation procedure.  
sim.tbs <- function(gen=list(n=100,cens=0,lambda=2,xi=3,beta=1,dist="norm"),prefix="",
                    n.copies=1,initial.seed=1234,method="Rsolnp",verbose=TRUE) {
  # gen: is a list with the parameters to generate from a TBS model
  #          n: sample size
  #       cens: censor rate
  #     lambda, xi, beta: TBS parameters
  # n.copies: how many copies will be performed in the simulation
  # initial.seed: to control the simulation

  result.1  <- matrix(gen$dist,n.copies,1)
  result.2  <- matrix(NA,n.copies,19)

  for (i in 1:n.copies) {
    # Generating data
    initial.time <- .gettime()
    seed <- initial.seed+i-1
    set.seed(seed)
    d <- NULL
    d$time  <- rtbs(gen$n,gen$lambda,gen$xi,gen$beta,dist=gen$dist)
    d$delta <- rep(1,gen$n)
    censor  <- quantile(d$time,probs=(1-gen$cens))
    d$delta <- ifelse(d$time > censor, 0, 1)
    d$time  <- ifelse(d$time > censor, censor, d$time)

    # Original Survival function
    orig.y  <- ptbs(d$time[d$delta == 1],gen$n,gen$lambda,gen$xi,gen$beta,dist=gen$dist)

    # Tbs fit
    tbs.fit  <- tbs.survreg.mle(Surv(d$time,d$delta) ~ 1,dist=gen$dist,method=method,verbose=verbose)

    # Result matrix
    if (tbs.fit$convergence) {
      tbs.y <- ptbs(d$time[d$delta == 1],lambda=tbs.fit$par[1],xi=tbs.fit$par[2],
                    beta=tbs.fit$par[3],dist=gen$dist)
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$lambda,gen$xi,gen$beta,tbs.fit$par,
                        (gen$lambda-tbs.fit$par[1])^2,(gen$xi-tbs.fit$par[2])^2,
                        (gen$beta-tbs.fit$par[3])^2,
                        max(abs(orig.y-tbs.y)),mean(abs(orig.y-tbs.y)),median(d$time),
                        exp(gen$beta),exp(tbs.fit$par[3]),(exp(tbs.fit$par[3])-exp(gen$beta))^2,
                        run.time)
    } else {
      run.time <- .gettime() - initial.time
      result.2[i,] <- c(seed,gen$n,gen$cens,gen$lambda,gen$xi,gen$beta,rep(NA,8),
                        median(d$time),rep(NA,3),run.time)
    }
    rm(d)
  }
  
  colnames(result.1) <- "dist"
  colnames(result.2) <- c("gen.seed","n","censor","orig.lambda","orig.xi","orig.beta",
                          "tbs.lambda","tbs.xi","tbs.beta","lambda.sqe","xi.sqe","beta.sqe",
                          "tbs.Max.AE","tbs.Mean.AE","orig.median","data.median",
                          "tbs.median","tbs.sqe","run.time")
  out <- data.frame(result.1,result.2)

  name <- paste(prefix,"Sim_TBS-",gen$dist,"_",gen$n,"_",gen$cens,"_",gen$lambda,
                "_",gen$xi,"_",gen$beta,".csv",sep="")
  write.csv(out,file=name)
#  return(out)
}






