

.gettime <- function() {
  t <- proc.time()
  ## note that t[1]+t[2] are the process time, not the real computer time,
  ## so we still take proper values even if the machine is running other stuff
  out <- (t[1]+t[2])/60
  names(out) <- NULL
  return(out)
}

psemi <- function(x,eta,peta) {
  if ((all(eta > 0)) && (all(peta >= 0)) && (all(peta <= 1))) {
    out <- c(crossprod(sapply(x,punif,min=-eta,max=eta),peta))
  }
  return(out)
}

dsemi <- function(x,eta,peta) {
  if ((all(eta > 0)) && (all(peta >= 0)) && (all(peta <= 1))) {
    out <- c(crossprod(sapply(x,dunif,min=-eta,max=eta),peta))
  }
  return(out)
}

g.lambda <- function(x,lambda) {
  return(sign(x)*(abs(x)^lambda)/lambda)
}

lik.semi <- function(par,time,delta,x=NULL,neta) {
  if (!is.null(x)) {
    if (is.matrix(x)) {
      nbeta <- length(x[1,])
      x1 <- x[delta == 1,]
      x0 <- x[delta == 0,]
    }
    else {
      nbeta <- length(x)
      x1 <- x[delta == 1]
      x0 <- x[delta == 0]
    }
  } else {
    nbeta <- 1
    x1 <- 1
    x0 <- 1
  }

  lambda <- par[1]
  beta   <- par[2:(nbeta+1)]
  eta    <- par[(nbeta+2):(nbeta+neta+1)]
  peta   <- par[(nbeta+neta+2):length(par)]

  out <- log(0)
  if ((all(time > 0)) && (lambda > 0) && (all(eta > 0)) &&
      (all(peta >= 0)) && (all(peta <= 1)) && (sum(peta) == 1)) {
    out <- 0
    t1 <- time[delta == 1]
    t0 <- time[delta == 0]
    if (!is.null(t1)) {
      if (is.matrix(x)) {
        e1 <- g.lambda(log(t1),lambda)-g.lambda(c(crossprod(t(x1),beta)),lambda)
      } else {
        e1 <- g.lambda(log(t1),lambda)-g.lambda(x1*beta,lambda)
      }
      out <- out + sum(log((1/t1)*(abs(log(t1))^(lambda-1))*dsemi(e1,eta,peta)))
    }
    if (!is.null(t0)) {
      if (is.matrix(x)) {
        e0 <- g.lambda(log(t0),lambda)-g.lambda(c(crossprod(t(x0),beta)),lambda)
      } else {
        e0 <- g.lambda(log(t0),lambda)-g.lambda(x0*beta,lambda)
      }
      out <- out + sum(log(psemi(e0,eta,peta)))
    }
  }
  return(out)
}

logpost.semi <- function(par,time,delta,x=NULL,neta) {
  if (!is.null(x)) {
    if (is.matrix(x)) {
      nbeta <- length(x[1,])
    }
    else {
      nbeta <- length(x)
    }
  } else {
    nbeta <- 1
  }

  lambda <- par[1]
  beta   <- par[2:(nbeta+1)]
  eta    <- par[(nbeta+2):(nbeta+neta+1)]
  V      <- par[(nbeta+neta+2):length(par)]

  out <- log(0)
  if ((all(time > 0)) && (lambda > 0) && (all(eta > 0)) &&
      (all(V >= 0)) && (all(V <= 1))) {
    out <- 0
    out <- out + dunif(lambda,min=0,max=4,log=TRUE)
    out <- out + sum(dnorm(beta,0,5,log=TRUE))
    out <- out + sum(dgamma((eta^2),2,1/2,log=TRUE))
    out <- out + sum(dbeta(V,1,1,log=TRUE))
    
    V[length(V)+1] <- 1
    peta <- rep(0,length(V))
    peta[1] <- V[1]
    for (i in 2:length(V)) {
      peta[i] <- V[i]*exp(sum(log(1-V[1:(i-1)])))
    }

    out <- out + lik.semi(par=c(lambda,beta,eta,peta),time=time,delta=delta,x=x,neta=neta)
  }

  if (is.nan(out) || is.na(out))
    out <- log(0)
  return(out)
}

#logpost.semi(c(1,1,1,2,3,0.3,0.3),rgamma(10,2,2),rep(1,10),1,3)


be.semi <- function(formula,neta=10,max.time=-1,
                    guess.beta,guess.lambda,
                    burn=1000,jump=2,size=1000,scale=1,
                    seed=1234) {
  initial.time <- .gettime()
  set.seed(seed)
  if(max.time <= 0) {
    ## user didn't define a timeout, so we set to a large number
    max.time <- 1e10
  }

  ## check the class of formula
  if (attributes(formula)$class != "formula")
    stop("A formula argument is required")

  ## record the call arguments
  Call  <- match.call()
  ## read the information from within the formula to populate the required variables
  mf <- model.frame(formula=formula)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  time <- y[,1]
  delta <- y[,2]
  x.k   <- dim(x)[2]
  n     <- dim(x)[1]
  ## check if delta is an indicator function
  if (any((delta != 0) & (delta != 1)))  {
    stop("It is only accepted uncensored or right censored data")
  }

  if (missing(guess.lambda)) {
    guess.lambda <- 1
  }
  if (missing(guess.beta)) {
    guess.beta <- rep(0,x.k)
  }
  guess.eta  <- seq(0.01,max(time)+1,(max(time)+1-0.01)/(neta-1))
  guess.V    <- rep(1/neta,(neta-1))

  out <- NULL
  out$call <- Call
  out$time <- time[delta == 1]

  ## perform a series of verifications for the given arguments of the function
  if (length(guess.lambda) != 1)
    stop("guess.lambda is not a scalar")
  if (guess.lambda <= 0)
    stop("guess.lambda must be a positive number")
  if (length(guess.beta) != x.k)
    stop("guess.beta length is not in accordance with the model specification")
  guess <- c(guess.lambda,guess.beta,guess.eta,guess.V)

  if ((!is.integer(burn)) && (burn <= 0)) 
    stop("burn must be a positive integer number")
  if ((!is.integer(jump)) && (jump <= 0)) 
    stop("jump must be a positive integer number")
  if ((!is.integer(size)) && (size <= 0)) 
    stop("size must be a positive integer number")
  if (scale <= 0)
    stop("scale must be a positive number")

  if (!is.null(x)) {
    if (is.matrix(x)) {
      if (length(time) != length(x[,1]))
        stop("length of time is different of length of x")
    } else {
      if (length(time) != length(x))
        stop("length of time is different of length of x")
      x <- matrix(x,length(x),1)
    }
    if(length(beta) > 1)
      beta <- matrix(beta,length(beta),1)
  } else {
    x <- matrix(1,length(time),1)
  }

  ## call the Metropolis algorithm for MCMC
  chain <- metrop(obj=logpost.semi,initial=guess,time=time,
                                      delta=delta,x=x,neta=neta,
                                      nbatch=(size-1)*jump+burn,blen=1,nspac=1,scale=scale)

  chain <- try(evalWithTimeout(metrop(obj=logpost.semi,initial=guess,time=time,
                                      delta=delta,x=x,neta=neta,
                                      nbatch=(size-1)*jump+burn,blen=1,nspac=1,scale=scale),
                               timeout=max.time*60,onTimeout="error"),silent=TRUE)
  if(length(class(chain)) == 1 && class(chain) == "try-error") {
    stop("Time limit exceeded")
  }
  out$post <- chain$batch[seq(burn,length(chain$batch[,1]),jump),]
  out$accept <- chain$accept

  # evaluating the point estimates
  out$run.time <- .gettime() - initial.time

  class(out) <- "tbs.survreg.be.semi"
  return(out)
}

#library(survival)
#data(colon)
#s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ 1,dist="norm",method=c("BFGS","Nelder-Mead"),verbose=TRUE)
b=be.semi(Surv(colon$time,colon$status==1) ~ 1,burn=10,jump=2,size=10,scale=0.05)




