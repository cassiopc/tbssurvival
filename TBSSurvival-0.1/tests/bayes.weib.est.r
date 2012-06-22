# Part of library 'TBSSurvival' for R (http://www.R-project.org)
#
# Copyright (C) 2010-2012, Adriano Polpo and D. Sinha.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
# 
# You should have received a copy of the GNU Library General Public
# License along with this library; A copy of the GNU General Public
# License is available on the World Wide Web at
# http://www.gnu.org/copyleft/gpl.html. You can also obtain it by
# writing to the Free Software Foundation, Inc., 51 Franklin St,
# Fifth Floor, Boston, MA  02110-1301, USA.

bayes.weib.est <- function(kick,time,delta,burn,jump,size,scale) {
  require("mcmc")
  require("coda")
  require("survival")
  initial.time <- .gettime()
  out <- NULL
  chain <- metrop(obj=lik.weib,initial=kick,time=time,delta=delta,nbatch=(size-1)*jump+burn,blen=1,nspac=1,scale=scale)
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

