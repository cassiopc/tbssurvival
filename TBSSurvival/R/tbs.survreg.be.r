# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
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

## TBS estimation using a Bayesian approach. The lack of a closed form
## solution forces us to tackle the problem by MCMC.
## max.time is a time limit in minutes (<= 0 means no limit)
## formula is specification containing a Surv model with right-censored data as in the package survival.
## dist defines the error distribution: "norm", "doubexp", "t", "cauchy", "logistic".
## guess.beta, guess.lambda, guess.xi are initial value of the Markov Chain (beta has to have the same
## of elements as covariates, lambda and xi are scalars.
## burn-in is the number of firsts samples of posterior to not use, and jump is the number of jump
## between each sample of posterior to avoid the problem of auto-correlation between the samples.
## size is the final sample size of the posterior. scale is a parameter of the `metrop' function which
## controls the acceptance rate. prior.mean and prior.sd define the parameters of the normal prior for
## the MCMC (by default, they are equal to 5 and 5).
## accept: fraction of Metropolis proposals accepted.
tbs.survreg.be <- function(formula,dist="norm",max.time=-1,
                           guess.beta,guess.lambda,guess.xi,
                           burn=1000,jump=2,size=1000,scale=1,
                           prior.mean=NULL,prior.sd=NULL,
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
  if (missing(guess.xi)) {
    guess.xi <- 1
  }
  if (missing(guess.beta)) {
    guess.beta <- rep(0,x.k)
  }

  out <- NULL
  out$call <- Call
#  if (is.matrix(x))
#    out$x <- x[order(time),]
#  else
#    out$x <- x[order(time)]
#  out$delta <- delta[order(time)]
#  out$time <- time[order(time)]
  out$time <- time[delta == 1]

  ## perform a series of verifications for the given arguments of the function
  if (length(guess.lambda) != 1)
    stop("guess.lambda is not a scalar")
  if (guess.lambda <= 0)
    stop("guess.lambda must be a positive number")
  if (length(guess.xi) != 1)
    stop("guess.xi is not a scalar")
  if (guess.xi <= 0)
    stop("guess.xi must be a positive number")
  if (length(guess.beta) != x.k)
    stop("guess.beta length is not in accordance with the model specification")
  guess <- c(guess.lambda,guess.xi,guess.beta)

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

  if (is.null(prior.mean)) {
    prior.mean <- 5
  } else {
    if (!is.vector(prior.mean))
      stop("mean is not a vector/scalar")
    if ((length(prior.mean) != 1) || (length(prior.mean) != length(guess[3:length(guess)]))) {
      stop(paste("length mean is different of 1 or ",length(guess[3:length(guess)]),sep=""))
    }
  }
  if (is.null(prior.sd)) {
    prior.sd <- 5
  } else {
    if (!is.vector(prior.sd))
      stop("sd is not a vector/scalar")
    if ((length(prior.sd) != 1) || (length(prior.sd) != length(guess[3:length(guess)]))) {
      stop(paste("length sd is different of 1 or ",length(par[3:length(guess)]),sep=""))
    }
    if (prior.sd <= 0)
      stop("prior.sd must be a positive number")
  }

  ## call the Metropolis algorithm for MCMC
  chain <- try(evalWithTimeout(metrop(obj=.logpost,initial=guess,time=time,delta=delta,dist=dist,x=x,
                                      mean=prior.mean,sd=prior.sd,
                                      nbatch=(size-1)*jump+burn,blen=1,nspac=1,scale=scale),
                               timeout=max.time*60,onTimeout="error"),silent=TRUE)
  if(length(class(chain)) == 1 && class(chain) == "try-error") {
    stop("Time limit exceeded")
  }
  out$post <- chain$batch[seq(burn,length(chain$batch[,1]),jump),]
  out$accept <- chain$accept

  # evaluating the point estimates
  if (x.k != 1) {
    par    <- c(mean(out$post[,1]),mean(out$post[,2]),apply(out$post[,3:length(out$post[1,])],2,mean))
    par.sd <- c(sd(out$post[,1]),    sd(out$post[,2]),apply(out$post[,3:length(out$post[1,])],2,sd))
  } else { 
    par    <- c(mean(out$post[,1]),mean(out$post[,2]),mean(out$post[,3:length(out$post[1,])]))
    par.sd <- c(sd(out$post[,1]),    sd(out$post[,2]),  sd(out$post[,3:length(out$post[1,])]))
  }
  out$lambda <- par[1]
  out$xi     <- par[2]
  out$beta   <- par[3:length(par)]
  out$lambda.sd <- par.sd[1]
  out$xi.sd     <- par.sd[2]
  out$beta.sd   <- par.sd[3:length(par.sd)]

  text <- c("lambda","xi")
  if (length(out$beta) == 1) {
    text <- c(text,"beta")
  } else {
    for (i in 1:length(out$beta))
      text <- c(text,paste("beta",i-1,sep=""))
  }
  colnames(out$post) <- text
  # evaluating the interval estimates
  par.HPD   <- cbind(c(HPDinterval(as.mcmc(out$post[,1]),0.95)),
                     c(HPDinterval(as.mcmc(out$post[,2]),0.95)))
  for (i in 3:length(out$post[1,])) {
    par.HPD <- cbind(par.HPD,
                     c(HPDinterval(as.mcmc(out$post[,i]),0.95)))
  }
  out$lambda.HPD <- par.HPD[,1]
  out$xi.HPD     <- par.HPD[,2]
  out$beta.HPD   <- par.HPD[,3:length(par.HPD[1,])]

  # evaluating DIC
  aux.loglik  <- rep(0,size)
  for (j in 1:size)
    aux.loglik[j] <- c(.lik.tbs(out$post[j,],time,delta,dist,x))
  loglik <- mean(-2*aux.loglik)
  rm(aux.loglik)
  out$DIC <- 2*loglik+2*.lik.tbs(par,time,delta,dist,x)

  # evaluating error of the model
  aux.error    <- matrix(0,length(time),size)
  for (j in 1:size) {
      if (x.k != 1) {
        aux.error[,j]    <- (.g.lambda(log(time),out$post[j,1])-
           .g.lambda(x%*%matrix(out$post[j,3:length(out$post[1,])],length(out$post[j,3:length(out$post[1,])]),1),out$post[j,1]))
      } else {
        aux.error[,j]    <- (.g.lambda(log(time),out$post[j,1])-
                             .g.lambda(out$post[j,3]*x,out$post[j,1]))
      }
  }
  error <- rep(0,length(time))
  for (i in 1:length(time))
    error[i] <- mean(aux.error[i,])
  rm(aux.error)
  out$error <- error[delta == 1]
  out$error.dist <- dist

  aux <- .test.tbs(out$lambda,out$xi,out$beta,x,dist,time,type="d")
  ## for the plot
  if (length(out$beta) == 1) {
    if (unique(aux$x[,1]) == 1) {
      out$x <- 1
      attr(out$x,"plot") <- 1
    } else if (length(unique(aux$x[,1]) <= 4)) {
      out$x <- unique(aux$x[,1])
      attr(out$x,"plot") <- 2
    }
  } else if ((length(out$beta) == 2) && (unique(aux$x[,1]) == 1) &&
             (length(unique(aux$x[,2])) <= 4)) {
    out$x <- unique(aux$x[,2])
    attr(out$x,"plot") <- 3
  }
  # time spent to do BE
  out$run.time <- .gettime() - initial.time

  class(out) <- "tbs.survreg.be"
  return(out)
}

print.tbs.survreg.be <- function(x, ...) {
  print(x$call)
  if (x$error.dist == "norm")
    text.dist <- "normal"
  if (x$error.dist == "t")
    text.dist <- "t-student"
  if (x$error.dist == "cauchy")
    text.dist <- "Cauchy"
  if (x$error.dist == "doubexp")
    text.dist <- "Double exponential"
  if (x$error.dist == "logistic")
    text.dist <- "logistic"
  cat("\n",sep="")
  cat("--------------------------------------------------------\n",sep="")
  cat(" TBS model with ",text.dist," error distribution (BE).\n",sep="")
  cat("\n",sep="")

  aux1 <- nchar(sprintf("%.4f",c(x$lambda,x$xi,x$beta)))
  aux2 <- nchar(sprintf("%.4f",c(x$lambda.sd,x$xi.sd,x$beta.sd)))
  text <- c("Estimates","Std. Error")
  auxc <- nchar(text)
  auxc1 <- c(ifelse(max(aux1) > auxc[1],abs(max(aux1)-auxc[1]),0),
             ifelse(max(aux2) > auxc[2],abs(max(aux2)-auxc[2]),0))
  auxc <- auxc+auxc1

  cat("         ",rep(" ",auxc1[1]+1),text[1],
                  rep(" ",auxc1[2]+1),text[2],"\n",sep="")
  cat("  lambda: ",
      rep(ifelse(auxc[1] > aux1[1]," ",""),abs(auxc[1]-aux1[1])),sprintf("%.4f",x$lambda)," ",
      rep(ifelse(auxc[2] > aux2[1]," ",""),abs(auxc[2]-aux2[1])),sprintf("%.4f",x$lambda.sd),
      "\n",sep="")
  cat("      xi: ",
      rep(ifelse(auxc[1] > aux1[2]," ",""),abs(auxc[1]-aux1[2])),sprintf("%.4f",x$xi)," ",
      rep(ifelse(auxc[2] > aux2[2]," ",""),abs(auxc[2]-aux2[2])),sprintf("%.4f",x$xi.sd),
      "\n",sep="")
  if (length(x$beta) == 1) {
    cat("    beta: ",
        rep(ifelse(auxc[1] > aux1[3]," ",""),abs(auxc[1]-aux1[3])),sprintf("%.4f",x$beta)," ",
        rep(ifelse(auxc[2] > aux2[3]," ",""),abs(auxc[2]-aux2[3])),sprintf("%.4f",x$beta.sd)," ",
        "\n",sep="")
  } else {
    for (i in 1:length(x$beta)) {
      if (i-1 < 10) {
        cat("   beta",i-1,": ",
            rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
            rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.sd[i])," ",
            "\n",sep="")
      } else {
        cat("  beta",i-1,": ",
            rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
            rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.sd[i])," ",
            "\n",sep="")
      }
    }
  }
  cat("\n",sep="")
  cat("     DIC: ",sprintf("%.4f",x$DIC),"\n",sep="")
  cat("\n",sep="")
  cat("Run time: ", sprintf("%.2f",x$run.time)," min \n",sep="")
  cat("--------------------------------------------------------\n",sep="")
  cat("\n",sep="")
}

summary.tbs.survreg.be <- function(x, ...) {
  if (x$error.dist == "norm")
    text.dist <- "normal"
  if (x$error.dist == "t")
    text.dist <- "t-student"
  if (x$error.dist == "cauchy")
    text.dist <- "Cauchy"
  if (x$error.dist == "doubexp")
    text.dist <- "Double exponential"
  if (x$error.dist == "logistic")
    text.dist <- "logistic"
  cat("--------------------------------------------------------\n",sep="")
  cat(" TBS model with ",text.dist," error distribution (BE).\n",sep="")
  cat("\n",sep="")

  auxt1 <- sprintf("%.4f",c(x$lambda,x$xi,x$beta))
  aux1  <- nchar(auxt1)
  auxt2 <- sprintf("%.4f",c(x$lambda.sd,x$xi.sd,x$beta.sd))
  aux2  <- nchar(auxt2)
  if (length(x$beta)) {
    auxt3 <- sprintf("%.4f",c(x$lambda.HPD[1],x$xi.HPD[1],x$beta.HPD[1]))
    aux3  <- nchar(auxt3)
    auxt4 <- sprintf("%.4f",c(x$lambda.HPD[2],x$xi.HPD[2],x$beta.HPD[2]))
    aux4  <- nchar(auxt4)
  } else {
    auxt3 <- sprintf("%.4f",c(x$lambda.HPD[1],x$xi.HPD[1],x$beta.HPD[1,]))
    aux3  <- nchar(auxt3)
    auxt4 <- sprintf("%.4f",c(x$lambda.HPD[2],x$xi.HPD[2],x$beta.HPD[2,]))
    aux4  <- nchar(auxt4)
  }

  text <- c("Estimates","Std. Error","95% HPD CI")
  text2 <- c("  lambda: ","      xi: ")
  if (length(x$beta) == 1) {
    text2 <- c(text2,"    beta: ")
  } else {
    for (i in 1:length(x$beta)) {
      if (i-1 < 10) {
        text2 <- c(text2,paste("   beta",i-1,": ",sep=""))
      } else {
        text2 <- c(text2,paste("  beta",i-1,": ",sep=""))
      }
    }
  }
  auxc <- nchar(text)

  auxc2    <- max(aux1,auxc[1]) -auxc[1]
  auxc2[2] <- max(aux2,auxc[2]) -auxc[2]
  auxc2[3] <- max(aux3+aux4+3,auxc[1]) -auxc[3]
  auxc2 <- ifelse(auxc2 < 0,0,auxc2)
  auxc[3] <- max(aux3)
  auxc[4] <- max(aux4)

  cat("         ",rep(" ",auxc2[1]+1),text[1],
                  rep(" ",auxc2[2]+1),text[2],
                  rep(" ",auxc2[3]+1),text[3],"\n",sep="")
  for (i in 1:length(text2)) {
    cat(text2[i],
        rep(ifelse(auxc[1] > aux1[i]," ",""),abs(auxc[1]-aux1[i])),auxt1[i]," ",
        rep(ifelse(auxc[2] > aux2[i]," ",""),abs(auxc[2]-aux2[i])),auxt2[i]," (",
        rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),auxt3[i],", ",
        rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),auxt4[i],")",
        "\n",sep="")
  }
  cat("\n",sep="")
  cat("Summary statistic of the posterior mean of the error\n")
  cat("for the TBS model:\n",sep="")
  print(summary(x$error))
  cat("\n",sep="")
  if (length(x$beta) == 1) {
    cat("Estimated quantiles of time event:\n",sep="")
    aux1 <- qtbs(c(0.05,0.25,0.5,0.75,0.95),
                 lambda=x$lambda,
                 xi=x$xi,
                 beta=x$beta,
                 dist=x$error.dist)
    aux2 <- nchar(trunc(aux1))
    cat("   ",rep(" ",aux2[1]),"5%   ",rep(" ",aux2[2]),"25%   ",rep(" ",aux2[3]),
        "50%   ",rep(" ",aux2[4]),"75%   ",rep(" ",aux2[5]),"95%\n",sep="")
    cat(sprintf("%.4f",aux1[1])," ",sprintf("%.4f",aux1[2])," ",sprintf("%.4f",aux1[3])," ",
        sprintf("%.4f",aux1[4])," ",sprintf("%.4f",aux1[5]),"\n",sep="")
  }  
  cat("--------------------------------------------------------\n",sep="")
  cat("\n",sep="")
}

plot.tbs.survreg.be <- function(x, HPD=TRUE, ...) {
  h <- 1000
  if (!exists("xlim")) {
    LB <- 0.0001
    UB <- max(x$time)*1.1
  } else {
    LB <- xlim[1]
    UB <- xlim[2]
  }
  if (!exists("xlab")) {
    xlab <- c("time","time","error")
  } else if (length(xlab) != 3) {
      stop("The 'xlab' length must be equal to 3")
  }
  if (!exists("ylab")) {
    ylab <- c("S(t)","h(t)","Density")
  } else if (length(ylab) != 3) {
      stop("The 'ylab' length must be equal to 3")
  }
  if (!exists("main")) {
    main <- c("Survival function (BE)","Hazard function (BE)","Posterior Mean of Error Distribution (BE)")
  } else if (length(main) != 3) {
      stop("The 'main' length must be equal to 3")
  }
  axis.t <- seq(LB,UB,(UB-LB)/(h-1))

  par(ask=TRUE)
  for (k in 1:2) {
    if (attr(x$x,"plot") == 1) {
      axis.y <- matrix(NA,length(axis.t),length(x$post[,1]))
      if (k == 1) { ### Survival plot
        for (j in 1:length(x$post[,1])) {
          axis.y[,j] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                               beta=x$post[j,3],dist=x$error.dist)
        }
        plot(axis.t,apply(axis.y,1,mean),type="l",xlim=c(0,UB),ylim=c(0,1),
             ylab=ylab[k],xlab=xlab[k],main=main[k], ...)
      } else { ### Hazard plot
        for (j in 1:length(x$post[,1])) {
          axis.y[,j] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                             beta=x$post[j,3],dist=x$error.dist)
        }
        plot(axis.t,apply(axis.y,1,mean),type="l",xlim=c(0,UB),ylab=ylab[k],
             xlab=xlab[k],main=main[k], ...)
      }
      if (HPD) {
        lines(axis.t,apply(axis.y,1,quantile,probs=0.025),type="l",col="gray50", ...)
        lines(axis.t,apply(axis.y,1,quantile,probs=0.975),type="l",col="gray50", ...)
      }
    } else if ((attr(x$x,"plot") == 2) || (attr(x$x,"plot") == 3)) {
      if (!exists("lty")) {
        lty <- seq(1,length(x$x),1)
      } else {
        if (length(lty) != length(x$x)) {
          stop(paste("The 'lty' length must be equal to ",length(x$x),sep=""))
        }
      }
      if (!is.numeric(col)) {
        col <- rep(1,length(x$x))
      } else {
        if (length(col) != length(x$x)) {
          stop(paste("The 'col' length must be equal to ",length(x$x),sep=""))
        }
      }
      if (!exists("lwd")) {
        lwd <- rep(1,length(x$x))
      } else {
        if (length(lwd) != length(x$x)) {
          stop(paste("The 'lwd' length must be equal to ",length(x$x),sep=""))
        }
      }
      axis.y <- array(NA,c(length(axis.t),length(x$post[,1]),length(x$x)))
      for (i in 1:length(x$x)) {
        if (attr(x$x,"plot") == 2) {
          if (k == 1) { ### Survival plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                   beta=x$post[j,3]*x$x[i],dist=x$error.dist)
            }
          } else { ### Hazard plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                 beta=x$post[j,3]*x$x[i],dist=x$error.dist)
            }
          }
        } else {
          if (k == 1) { ### Survival plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                     beta=(x$post[j,3]+x$post[j,4]*x$x[i]),dist=x$error.dist)
            }
          } else { ### Hazard plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                   beta=(x$post[j,3]+x$post[j,4]*x$x[i]),dist=x$error.dist)
            }
          }
        }
        if (i == 1) {
          if (k == 1) { ### Survival plot
            plot(axis.t,apply(axis.y[,,i],1,mean),xlim=c(0,UB),ylim=c(0,1),ylab=ylab[k],xlab=xlab[k],main=main[k],
                 type="l",lty=lty[i],lwd=lwd[i],col=col[i], ...)
          } else {
            plot(axis.t,apply(axis.y[,,i],1,mean),xlim=c(0,UB),ylab=ylab[k],xlab=xlab[k],main=main[k],
                 type="l",lty=lty[i],lwd=lwd[i],col=col[i], ...)
          }
        } else {
          lines(axis.t,apply(axis.y[,,i],1,mean),lty=lty[i],lwd=lwd[i],col=col[i])
        }
        if (HPD) {
          lines(axis.t,apply(axis.y[,,i],1,quantile,probs=0.025),type="l",col="gray50", ...)
          lines(axis.t,apply(axis.y[,,i],1,quantile,probs=0.975),type="l",col="gray50", ...)
        }
      }
    }
  }

  ## Histogram for error distribution:
  bound <- max(abs(c(min(x$error),max(x$error))))
  hist(x$error,probability=TRUE,xlim=c(-bound,bound),main=main[3],xlab=xlab[3],ylab=ylab[3])
  axis.x <- seq(-bound,bound,2*bound/(h-1))
  axis.y <- matrix(NA,length(axis.x),length(x$post[,1]))
  for (j in 1:length(x$post[,1])) {
    axis.y[,j] <- .choice(axis.x,xi=x$post[j,2],dist=x$error.dist,type="d")
  }
  lines(axis.x,apply(axis.y,1,mean))
  if (HPD) {
    lines(axis.x,apply(axis.y,1,quantile,probs=0.025),type="l",col="gray50", ...)
    lines(axis.x,apply(axis.y,1,quantile,probs=0.975),type="l",col="gray50", ...)
  }

  if (FALSE) {
    ## Q-Q plot for error distribution:
    axis.x <- matrix(NA,length(x$error),length(x$post[,1]))
    for (j in 1:length(x$post[,1])) {
      axis.x[,j] <- qqplot(.choice(ppoints(length(x$error)), xi=x$post[j,2],
                                   dist=x$error.dist, type="q"),x$error,plot.it=FALSE)$x
    }
    axis.y <- qqplot(.choice(ppoints(length(x$error)), xi=x$xi,
                             dist=x$error.dist, type="q"),x$error,plot.it=FALSE)$y
    plot(apply(axis.x,1,mean),axis.y,
         main = expression("Q-Q plot for error"),xlab="Theoretical Quantiles",ylab="Sample Quantiles")
    qqline(x$error, distribution = function(p) {
                                     fn.aux <- function(p) {
                                       aux <- rep(NA,length(x$post[,1]))
                                       for (j in 1:length(x$post[,1]))
                                         aux[j] <- .choice(p, xi=x$post[j,2], dist=x$error.dist, type="q")
                                       return(mean(aux))
                                     }
                                     aux <- Vectorize(fn.aux,"p")
                                     return(aux(p))
                                   },
           prob = c(0.25, 0.75), col = 2, lwd=2)
  }

  #autocorrelation plot
  acf(x$post[,1],main=expression(lambda))
  acf(x$post[,2],main=expression(xi))
  if (length(x$post[1,]) == 3) {
    acf(x$post[,3],main=expression(beta))
  } else {
    for (i in 3:length(x$post[1,])) {
      acf(x$post[,i],main=bquote(beta[.(i-3)]))
    }
  }

  #time series plot
  plot(ts(x$post[,1]),ylab=expression(lambda),xlab="i")
  plot(ts(x$post[,2]),ylab=expression(xi),xlab="i")
  if (length(x$post[1,]) == 3) {
    plot(ts(x$post[,3]),ylab=expression(beta),xlab="i")
  } else {
    for (i in 3:length(x$post[1,])) {
      plot(ts(x$post[,i]),ylab=bquote(beta[.(i-3)]),xlab="i")
    }
  }

  par(ask=FALSE)
}








