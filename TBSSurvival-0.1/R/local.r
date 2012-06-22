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

#######################################################################
## This file contains
## Auxiliar functions - NOT to be used directly by user.
#######################################################################

.onAttach <- function(lib,pkg)
{
  packageStartupMessage("TBSSurvival 0.1 loaded\n")
}

##  Density for misture of uniforme-exponential. This density is use as prior for GPT model, has between (a,b)
##  an uniform distribution and after `b' the tail is exponential. the parameter `p' volume of uniform part and
##  `1-p' the volume of the tail.
##  \item{x}{vector of quantiles.}
##  \item{p}{number `0 < p < 1'.}
##  \item{a}{parameter of uniforme, default=0.}
##  \item{b}{parameter that difine the end of uniforme part.}
## \value{ `dunif.exp' gives the density. }
.dunif.exp <- function(x,a=0,b,p) {
  xb <- 1*(x <= b)
  out <- xb*p*dunif(x,a,b) + (1-xb)*dexp(x[i],-log(1-p)/b)

#  out <- rep(0,length(x))
#  for (i in 1:length(x)) {
#    if (x[i] <= b)
#      out[i] <- (x[i] <= b) p*dunif(x[i],a,b)
#    else
#      out[i] <- dexp(x[i],-log(1-p)/b)
#  }
  return(out)
}

.logpost <- function(par,time,delta,dist,x=NULL,mean=5,sd=5) {
  if ((par[1] > 0) && (par[2] > 0)) {
    if (dist != "t") {
      out <- log(.dunif.exp(par[1],a=0.00001,b=3,p=0.8))+
             log(.dunif.exp(par[2],a=0.00001,b=2,p=0.9))+
             sum(log(dnorm(par[3:length(par)],mean,sd)))+
             .lik.tbs(par=par,time=time,delta=delta,dist=dist,x=x)
    } else {
      out <- log(.dunif.exp(par[1],a=0.00001,b=3,p=0.8))+
             sum(log(dnorm(par[3:length(par)],mean,sd)))+
             .lik.tbs(par=par,time=time,delta=delta,dist=dist,x=x)
    }
  } else {
    out <- log(0)
  }
  
  if (is.nan(out) || is.na(out))
    out <- log(0)
  return(out)
}

.gettime <- function() {
  ## this function returns the current amount of spent time by
  ## the current R process (user time + system time) in minutes
  t <- proc.time()
  return((t[1]+t[2])/60)
}

.lik.tbs <- function(par,time,delta,dist="norm",x=NULL,notinf=FALSE)
{
  lambda <- par[1]
  xi     <- par[2]
  beta   <- par[3:length(par)]

  out <- log(0)
  if ((xi > 0) && (all(time > 0)) && (lambda > 0))
  {
    if (is.matrix(x)) {
      d.aux <- dtbs(time=time[delta==1],lambda=lambda,xi=xi,beta=beta,x=x[delta==1,],dist=dist)
      if (min(delta) == 0) {
        s.aux <- 1-ptbs(time=time[delta==0],lambda=lambda,xi=xi,beta=beta,x=x[delta==0,],dist=dist)
        out <- sum(log(d.aux))+sum(log(s.aux))
      } else {
        out <- sum(log(d.aux))
      }
    } else {
      d.aux <- dtbs(time=time[delta==1],lambda=lambda,xi=xi,beta=beta,x=x[delta==1],dist=dist)
      if (min(delta) == 0) {
        s.aux <- 1-ptbs(time=time[delta==0],lambda=lambda,xi=xi,beta=beta,x=x[delta==0],dist=dist)
        out <- sum(log(d.aux))+sum(log(s.aux))
      } else {
        out <- sum(log(d.aux))
      }
    }
  }
  if(notinf && out==-Inf) return(-10e10)
  return(out)
}

.tbs.survreg <- function(formula,dist="norm",method="BFGS",kick=NULL,nstart=4,verbose=FALSE) {
  initial.time <- .gettime()
  require("survival")

  if (attributes(formula)$class != "formula")
    stop("A formula argument is required")

  ## read the information from within the formula to populate the required variables
  mf <- model.frame(formula=formula)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  time <- y[,1]
  delta <- y[,2]
##  time  <- eval(attr(terms(formula),"variables"))[[1]][,1]
##  delta <- eval(attr(terms(formula),"variables"))[[1]][,2]
##  x     <- model.matrix(formula)
  x.k   <- dim(x)[2]
  n     <- dim(x)[1]
  if (any((delta != 0) & (delta != 1)))  {
    stop("Only uncesored or right censored data are allowed")
  }

  out <- NULL

  ## check if starting point was given or not, and build one in case
  if (is.null(kick)) {
    nparam <- 2
    if (!is.null(x)) {
      if (is.matrix(x))
        nparam <- nparam+length(x[1,])
      else
        nparam <- nparam+1
    }
    ## betas can be anything, we sample uniformly from -10 to 10
    kick <- 20*runif(nparam)-10
    ## lambda and xi have to be positive, and "desirable" values are not very high...
    kick[1] <- 5*runif(1)+0.0001 ## lambda
    kick[2] <- 10*runif(1)+0.0001 ## xi
  } else {
    nparam <- length(kick)
  }

  if(method=="Rsolnp") {
    if(library("Rsolnp",quietly=TRUE,logical.return=TRUE)==FALSE) {
      out$method <- "Rsolnp: not installed"
      out$convergence <- FALSE
      return(out)
    } else {
      ## Rsolnp needs some bounds for the unknowns. We arbitrarily define them to be between -100 and 100.
      LB = rep(-20,nparam)
      UB = rep(20,nparam)
      ## unless for lambda and xi, which are positive
      LB[1] = 0.0001
      LB[2] = 0.0001
      ## upper bound for xi is not very clear, but 1000 should be enough
      UB[2] = 1000
      ## try to run the solver, using parallel computing if available
      if(require("multicore",quietly=TRUE)) {
        if(verbose) cat('RSOLNP-multicore: ')
        ans = try(gosolnp(pars = NULL, fixed = NULL, fun = function(pars, n) { -.lik.tbs(pars,time=time,delta=delta,x=x,dist=dist,notinf=TRUE) },
          LB = LB, UB = UB, control = list(outer.iter = 100, trace = 0),
          distr = rep(1, length(LB)), distr.opt = list(), n.restarts = nstart, n.sim = 200, parallel=TRUE,parallel.control=list(pkg="multicore",core=4), rseed = runif(n=1,min=1,max=1000000), n = nparam))
      } else {
        if(require("snowfall",quietly=TRUE)) {
          if(verbose) cat('RSOLNP-snowfall: ')
          ans = try(gosolnp(pars = NULL, fixed = NULL, fun = function(pars, n) { -.lik.tbs(pars,time=time,delta=delta,x=x,dist=dist,notinf=TRUE) },
            LB = LB, UB = UB, control = list(outer.iter = 100, trace = 0),
            distr = rep(1, length(LB)), distr.opt = list(), n.restarts = nstart, n.sim = 200, parallel=TRUE,parallel.control=list(pkg="snowfall",core=4), rseed = runif(n=1,min=1,max=1000000), n = nparam))
        } else {
          if(verbose) cat('RSOLNP: ')
          ans = try(gosolnp(pars = NULL, fixed = NULL, fun = function(pars, n) { -.lik.tbs(pars,time=time,delta=delta,x=x,dist=dist,notinf=TRUE) },
            LB = LB, UB = UB, control = list(outer.iter = 100, trace = 0),
            distr = rep(1, length(LB)), distr.opt = list(), n.restarts = nstart, n.sim = 200, rseed = runif(n=1,min=1,max=1000000), n = nparam))
        }
      }
      if (class(ans) != "try-error") {
        ## process the solution in case one was found
        ## get parameters
        out$par <- ans$pars
        options(warn = -1)
        ## compute the std.error
        aux <- try(sqrt(diag(solve((ans$hessian)))),silent=TRUE)
        options("warn" = 0)
        out$std.error <- rep(NA,nparam)
        if (class(aux) != "try-error")
          out$std.error <- aux
        ## get the log-lik value
        out$log.lik <- -ans$values[length(ans$values)]
        if(verbose) cat(out$log.lik,'PARS:',ans$pars,'TIME:',ans$elapsed,'\n')
        out$error.dist <- dist
        ## compute error distances
        out$AIC  <- 2*nparam-2*out$log.lik
        out$AICc <- 2*nparam-2*out$log.lik + 2*nparam*(nparam+1)/(length(time)-nparam-1)
        out$BIC  <- -2*out$log.lik+nparam*log(length(time))
        out$method <- method
        out$convergence <- ans$convergence==0
        out$run.time <- .gettime() - initial.time
        return(out)
      } else {
        if(verbose) cat(' failed - fallback to BFGS\n')
        method="BFGS"
      }
    }
  }
  
  i <- 1
  est=NA
  ii=1
  if(verbose) cat(method,': ',sep='')
  inimethod=method
  repeat {
    valik=.lik.tbs(kick,time=time,delta=delta,x=x,dist=dist)
    if(!is.na(valik) && valik>-Inf) {
      aux <- try(optim(kick, .lik.tbs, time=time, delta=delta, dist=dist, x=x,
                       method=inimethod, control=list(fnscale=-1), hessian=TRUE),silent=TRUE)
      if (class(aux) != "try-error") {
        repeat {
          aux1 <- try(optim(aux$par, .lik.tbs, time=time, delta=delta, dist=dist, x=x,
                           method=method, control=list(fnscale=-1), hessian=TRUE),silent=TRUE)
          if (class(aux1) != "try-error") {
            if (aux1$value < aux$value + 0.0001) {
              ## 0.0001 is only for numerical reasons. Note that 0.0001 in the log value is anyway very very small...
              break
            }
            aux=aux1
          } else {
            break
          }
        }
      }
      if (class(aux) != "try-error") {
        if(is.na(est) || aux$value > est$value) {
          est = aux
          inimethod=method
        }
        i = i + 1
        if(verbose) cat('@')
      } else {
        if(verbose) cat('*')
        ii = ii + 1
      }
    } else {
      if(verbose) cat('.')
      ii = ii + 1
    }
    ## betas can be anything, we sample uniformly from -10 to 10
    kick <- 20*runif(nparam)-10
    ## lambda and xi have to be positive, and "desirable" values are not very high...
    kick[1] <- 5*runif(1)+0.0001 ## lambda
    kick[2] <- 10*runif(1)+0.0001 ## xi

    ## if enough starts have been tried, stop. Also stop if too many unsuccessfull tries have been made :(
    if(ii>100 && is.na(est)) {
      inimethod="SANN"
    }
    if(i>nstart || ii >500) {
      break
    }
    next
  }

  out$method <- method
  if(!is.na(est) && est$value > -Inf) {
    out$par <- est$par
    options(warn = -1)
    aux <- try(sqrt(diag(solve(-(est$hessian)))),silent=TRUE)
    options("warn" = 0)
    out$std.error <- rep(NA,nparam)
    out$run.time <- .gettime() - initial.time

    if (class(aux) != "try-error")
      out$std.error <- aux
    out$log.lik <- est$value
    out$error.dist <- dist
    if(verbose) cat(' ',out$log.lik,'PARS:',out$par,'TIME:',out$run.time,'\n')

    out$AIC  <- 2*nparam-2*est$value
    out$AICc <- 2*nparam-2*est$value + 2*nparam*(nparam+1)/(length(time)-nparam-1)
    out$BIC  <- -2*est$value+nparam*log(length(time))
    out$convergence <- TRUE
  } else {
    if(verbose) cat(' failed\n')
    out$run.time <- .gettime() - initial.time
    out$convergence <- FALSE
    warning("It was not possible to find a starting value")
  }
  return(out)
}

.test.tbs <- function(lambda, xi, beta, x=NULL, dist, time=NULL, type=NULL, p=NULL, n=NULL) {
  out   <- NULL
  out$x <- x
  out$beta <- beta
  if ((dist != "norm") && (dist != "doubexp") && (dist != "cauchy") &&
      (dist != "t")    && (dist != "logistic"))
    stop("dist: Distribution not available")
  if (!is.double(xi))
    stop("xi is not a number")
  if (is.matrix(xi))
    stop("xi is matrix")
  if ((!is.double(lambda)) || (length(lambda) != 1))
    stop("lambda is not a number or length != 1")
  if (!is.double(beta))
    stop("beta is not a (vector) number")
  if (is.matrix(beta))
    stop("beta is matrix")
  if (!is.null(x)) {
    if (is.matrix(x)) {
      if (length(beta) != length(x[1,]))
        stop(paste("size of beta != ",length(x[1,]),sep=""))
    }
    else {
      if ((length(beta) != 1) && (length(beta) != length(x)))
        stop("size of beta is not conform")
    }
  }
  else {
    if (length(beta) > 1)
      stop("x is wrong or length(beta) > 1")  
  }
  if (lambda <= 0)
    stop("lambda <= 0")
  if (xi <= 0)
    stop("xi <= 0")

  if (!is.null(type)) {
    if ((type == "d") || (type == "p")) {
      if (!is.double(time))
        stop("time is not a (vector) number")
      if (is.matrix(time))
        stop("time is matrix")
      if (any(time <= 0))
        stop("time <= 0")
      if (!is.null(x)) {
        if (is.matrix(x)) {
          if (length(time) != length(x[,1]))
            stop("length of time is different of length of x")
        }
        else {
          if (length(beta) == length(x)) {
            out$x <- matrix(x,1,length(x))
          } else {
            if (length(time) != length(x))
              stop("length of time is different of length of x")
            out$x <- matrix(x,length(x),1)
          }
        }
        out$beta <- matrix(beta,length(beta),1)
      }
      else {
        out$x <- matrix(1,length(time),1)
      }
    } else {
      if (type == "q") {
        if (!is.double(p))
          stop("p is not a (vector) number")
        if (is.matrix(p))
          stop("p is matrix")
        if (min(p) < 0)
          stop("p < 0")
        if (max(p) > 1)
          stop("p > 1")
      } else if (type == "r") {
          if (!is.double(n))
            stop("n is not a number")
          if (n %% 1 != 0)
            stop("n is not a integer number")
        }
      if (is.null(x)) {
        if (length(beta) > 1)
          stop("If x is omitted then beta must have length 1")
        out$x <- 1
      }
    }
  }

  return(out)
}

.choice <- function(x, xi, dist, type) {
  switch(dist,
         norm     = switch(type,
           d = dnorm(x,mean=0,sd=sqrt(xi)),
           p = pnorm(x,mean=0,sd=sqrt(xi)),
           q = qnorm(x,mean=0,sd=sqrt(xi)),
           r = rnorm(x,mean=0,sd=sqrt(xi))),
         t        = switch(type,
           d = dt(x,df=xi),
           p = pt(x,df=xi),
           q = qt(x,df=xi),
           r = rt(x,df=xi)),
         cauchy   = switch(type,
           d = dcauchy(x,location=0,scale=xi),
           p = pcauchy(x,location=0,scale=xi),
           q = qcauchy(x,location=0,scale=xi),
           r = rcauchy(x,location=0,scale=xi)),
         doubexp  = switch(type,
           d = dnormp(x,sigmap=xi,mu=0,p=1),
           p = pnormp(x,sigmap=xi,mu=0,p=1),
           q = qnormp(x,sigmap=xi,mu=0,p=1),
           r = rnormp(x,sigmap=xi,mu=0,p=1)),
         logistic = switch(type,
           d = dlogis(x,location=0,scale=xi),
           p = plogis(x,location=0,scale=xi),
           q = qlogis(x,location=0,scale=xi),
           r = rlogis(x,location=0,scale=xi)))
}


## \code{.g.lambda} gives the generalized power transformation function
.g.lambda <- function(x,lambda) {
  return(sign(x)*(abs(x)^lambda)/lambda)
}

## \code{.g.lambda.inv} the inverse of generalized power transformation function.
.g.lambda.inv <- function(x,lambda) {
  return(sign(x)*(abs(x*lambda)^(1/lambda)))
}
