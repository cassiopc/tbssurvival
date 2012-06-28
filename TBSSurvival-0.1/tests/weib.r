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

hweib <- function(x,shape,scale)
  return((shape/scale)*((x/scale)^(shape-1)))

dweib <- function(x,alpha,beta,mu=0)
  return((beta/alpha)*(((x-mu)/alpha)^(beta-1))*exp(-((x-mu)/alpha)^beta))

pweib <- function(x,alpha,beta,mu=0)
  return(1-exp(-((x-mu)/alpha)^beta))

qweib <- function(p,alpha,beta,mu=0)
  return(alpha*(-log(1-p))^(1/beta)+mu)

rweib <- function(n,alpha,beta,mu=0)
  return(alpha*(-log(runif(1)))^(1/beta)+mu)

lweib <- function(alpha,beta,mu=0,time,delta)
{
  shape <- beta
  scale <- alpha
  out <- log(0)
  if ((scale > 0) && (min(time) > mu) && (shape > 0))
  {
    d.aux <- shape*(scale^(-shape))*((time[delta == 1]-mu)^(shape-1))*exp(-(((time[delta == 1]-mu)/scale)^shape))
    out <- sum(log(d.aux))    
    if (length(time[delta == 0]) != 0) {
      s.aux <- exp(-(((time[delta == 0]-mu)/scale)^shape))
      out <- out + sum(log(s.aux))
    }
  }
  return(out)
}

bweib <- function(guess,time,delta,burn,jump,size,scale) {
  require("mcmc")
  require("coda")
  require("survival")
  initial.time <- .gettime()
  out <- NULL
  chain <- metrop(obj=lik.weib,initial=guess,time=time,delta=delta,nbatch=(size-1)*jump+burn,blen=1,nspac=1,scale=scale)
  out$post <- chain$batch[seq(burn,length(chain$batch[,1]),jump),]

  km <- survfit(formula = Surv(time, delta == 1) ~ 1)
  out$time <- km$time
  aux.hazard   <- matrix(0,length(out$time),size)
  aux.survival <- matrix(0,length(out$time),size)
  aux.density  <- matrix(0,length(out$time),size)
  for (j in 1:size) {
      aux.hazard[,j]   <- c(hazard.weibull(out$time,shape=out$post[j,1],scale=out$post[j,2]))
      aux.survival[,j] <-     c(1-pweibull(out$time,shape=out$post[j,1],scale=out$post[j,2]))
      aux.density[,j]  <-       c(dweibull(out$time,shape=out$post[j,1],scale=out$post[j,2]))
  }
  hazard   <- matrix(0,length(out$time),8)
  survival <- matrix(0,length(out$time),8)
  density  <- matrix(0,length(out$time),8)
  for (i in 1:length(out$time)) {
    hazard[i,]   <- c(  summary(aux.hazard[i,]),HPDinterval(as.mcmc(  aux.hazard[i,]),0.95))
    survival[i,] <- c(summary(aux.survival[i,]),HPDinterval(as.mcmc(aux.survival[i,]),0.95))
    density[i,]  <- c( summary(aux.density[i,]),HPDinterval(as.mcmc( aux.density[i,]),0.95))
  }
  rm(aux.hazard,aux.survival,aux.density)
  out$KS <- max(abs(survival[,3]-km$surv))

  aux.sum <- 0
  aux <- time[delta == 1]
  for (i in 1:length(aux))
    aux.sum <- aux.sum + log(density[out$time == aux[i],3])
  if (length(time[delta == 0]) != 0) {
    aux <- time[delta == 0]
    for (i in 1:length(aux))
      aux.sum <- aux.sum + log(survival[out$time == aux[i],3])
  }
  aux.sum <- -2*aux.sum
  rm(aux,density)

  aux.loglik  <- rep(0,size)
  for (j in 1:size)
      aux.loglik[j] <- c(lik.weib(out$post[j,],time,delta))
  loglik <- mean(-2*aux.loglik)
  rm(aux.loglik)

  out$DIC <- 2*loglik-aux.sum
  rm(aux.sum)

  out$par           <- c(median(out$post[,1]),median(out$post[,2]))
  out$par.std.error <- c(    sd(out$post[,1]),    sd(out$post[,2]))
  out$par.HPD       <- cbind(c(HPDinterval(as.mcmc(out$post[,1]),0.95)),
                             c(HPDinterval(as.mcmc(out$post[,2]),0.95)))
  out$par.KS <- max(abs((1-pweibull(out$time,shape=out$par[1],scale=out$par[2])) - km$surv))
  out$par.DIC <- 2*loglik+2*lik.weib(out$par,time,delta)

  out$survival <- survival
  out$hazard   <- hazard
  out$run.time <- .gettime() - initial.time
  return(out)
}

