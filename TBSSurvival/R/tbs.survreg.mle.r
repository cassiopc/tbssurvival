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

## Maximum likelihood estimation for TBS
## By default, try all optimization methods listed below and keep the best solution
## max.time is the time limit in minutes for each method (<= 0 means no limit), ncore the number of cores
##   to use in case multicore is available, nstart the number of feasible starting
##   points to use, dist has to be one of the available distributions, currently
##   one of "norm", "t", "cauchy", "doubexp", "logistic"
## formula is a R formula with a Surv object on the left side
tbs.survreg.mle <- function(formula,dist="all",method=c("BFGS","Rsolnp","Nelder-Mead","CG","SANN"),
                            verbose=FALSE,nstart=10,max.time=-1,ncore=1,seed=1234) {
  ## this meta-method only records the elapsed time and call the max likelihood estimation function
  ## for each of the methods given until one of them converges. It is supposed that at least one method
  ## is given, and that dist is one of those implemented by tbs.survreg.
  initial.time <- .gettime()

  set.seed(seed)
  fn.aux <- function(formula,dist,method,verbose,nstart,max.time,ncore) {
    initial.time <- .gettime()
    bestout <- NULL
    for(i in 1:length(method)) {
      ## call the estimation function
      out <- .tbs.survreg(formula,dist=dist,method=method[i],verbose=verbose,ncore=ncore,max.time=max.time,nstart=nstart)
      ## if converged, we are happy
      if(out$convergence) {
        if(is.null(bestout) || out$log.lik > bestout$log.lik) {
          bestout <- out
        }
      }
    }
    ## check if at least one method found a solution...
    if (is.null(bestout)) {
      if(verbose) cat('\nTBS: No method has found a solution\n')
      bestout$convergence <- FALSE
      bestout$method <- NULL
    }
    ## record the call arguments and the formula
    bestout$call <- match.call()
    bestout$formula <- formula
    ## run.time is returned in minutes
    bestout$run.time <- .gettime() - initial.time
  
    class(bestout) <- "tbs.survreg.mle"
    return(bestout)
  }

  if ((length(dist) > 1) || (dist == "all")) {
    out <- NULL
    if ((length(dist) == 1) && (dist == "all")) {
      dist <- c("norm","t","cauchy","doubexp","logistic")
    }
    aux <- Inf
    for (i in 1:length(dist)) {
      if (verbose) cat(dist[i],":\n",sep="")
      eval(parse(text=paste("out$",dist[i]," <- 1",sep=""))) 
      out[[i]] <- fn.aux(formula=formula,dist[i],method,verbose,nstart,max.time,ncore)
      if ((aux > out[[i]]$AIC) && (out[[i]]$convergence)) {
        aux <- out[[i]]$AIC
        best <- dist[i]
      }
    }
    if (exists("best")) {
      out$best <- best
    } else {
      out$best <- "none"
    }
  } else {
    out <- fn.aux(formula,dist,method,verbose,nstart,max.time,ncore)
  }
  
  out$run.time <- .gettime() - initial.time
  return(out)
}

print.tbs.survreg.mle <- function(x, ...) {
  print(x$call)
  if (x$convergence) {
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
    cat(" TBS model with ",text.dist," error distribution (MLE).\n",sep="")
    cat("\n",sep="")
  
    aux1 <- nchar(sprintf("%.4f",c(x$lambda,x$xi,x$beta)))
    aux2 <- nchar(sprintf("%.4f",c(x$lambda.std.error,x$xi.std.error,x$beta.std.error)))
    text <- c("Estimates","Std. Error")
    auxc <- nchar(text)
    auxc1 <- c(ifelse(max(aux1) > auxc[1],abs(max(aux1)-auxc[1]),0),
               ifelse(max(aux2) > auxc[2],abs(max(aux2)-auxc[2]),0))
    auxc <- auxc+auxc1
  
    cat("         ",rep(" ",auxc1[1]+1),text[1],
                    rep(" ",auxc1[2]+1),text[2],"\n",sep="")
    cat("  lambda: ",
        rep(ifelse(auxc[1] > aux1[1]," ",""),abs(auxc[1]-aux1[1])),sprintf("%.4f",x$lambda)," ",
        rep(ifelse(auxc[2] > aux2[1]," ",""),abs(auxc[2]-aux2[1])),sprintf("%.4f",x$lambda.std.error),
        "\n",sep="")
    cat("      xi: ",
        rep(ifelse(auxc[1] > aux1[2]," ",""),abs(auxc[1]-aux1[2])),sprintf("%.4f",x$xi)," ",
        rep(ifelse(auxc[2] > aux2[2]," ",""),abs(auxc[2]-aux2[2])),sprintf("%.4f",x$xi.std.error),
        "\n",sep="")
    if (length(x$beta) == 1) {
      cat("    beta: ",
          rep(ifelse(auxc[1] > aux1[3]," ",""),abs(auxc[1]-aux1[3])),sprintf("%.4f",x$beta)," ",
          rep(ifelse(auxc[2] > aux2[3]," ",""),abs(auxc[2]-aux2[3])),sprintf("%.4f",x$beta.std.error)," ",
          "\n",sep="")
    } else {
      for (i in 1:length(x$beta)) {
        if (i-1 < 10) {
          cat("   beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.std.error[i])," ",
              "\n",sep="")
        } else {
          cat("  beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.std.error[i])," ",
              "\n",sep="")
        }
      }
    }
    cat("\n",sep="")
    cat("     AIC: ",sprintf("%.4f",x$AIC),"\n",sep="")
    cat("    AICc: ",sprintf("%.4f",x$AICc),"\n",sep="")
    cat("     BIC: ",sprintf("%.4f",x$BIC),"\n",sep="")
    cat("\n",sep="")
    cat("Run time: ", sprintf("%.2f",x$run.time)," min \n",sep="")
    cat("--------------------------------------------------------\n",sep="")
    cat("\n",sep="")
  } else {
    cat("There was not obtained the convergence.\n",sep="")
  }
}

summary.tbs.survreg.mle <- function(x, ...) {
  if (x$convergence) {
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
    cat(" TBS model with ",text.dist," error distribution (MLE).\n",sep="")
    cat("\n",sep="")
  
    z.value <- rep(NA,length(x$beta))
    p.value <- rep(NA,length(x$beta))
    for (i in 1:length(x$beta)) {
      z.value[i] <- x$beta[i]/x$beta.std.error[i]
      p.value[i] <- 2*pnorm(-abs(z.value[i]),0,1)
    }
  
    aux1 <- nchar(sprintf("%.4f",c(x$lambda,x$xi,x$beta)))
    aux2 <- nchar(sprintf("%.4f",c(x$lambda.std.error,x$xi.std.error,x$beta.std.error)))
    aux3 <- nchar(sprintf("%.4f",z.value))
    p.value2 <- rep(NA,length(x$beta))
    for (i in 1:length(x$beta)) {
      if (p.value[i] < 0.0001) {
        p.value2[i] <- "< 0.0001"
      } else {
        p.value2[i] <- sprintf("%.4f",p.value[i])
      }
    }
    aux4 <- nchar(p.value2)
    text <- c("Estimates","Std. Error","z value","Pr(>|z|)")
    auxc <- nchar(text)
  
    auxc1 <- c(ifelse(max(aux1) > auxc[1],abs(max(aux1)-auxc[1]),0),
               ifelse(max(aux2) > auxc[2],abs(max(aux2)-auxc[2]),0),
               ifelse(max(aux3) > auxc[3],abs(max(aux3)-auxc[3]),0),
               ifelse(max(aux4) > auxc[4],abs(max(aux4)-auxc[4]),0))
  
    cat("         ",rep(" ",auxc1[1]+1),text[1],
                    rep(" ",auxc1[2]+1),text[2],
                    rep(" ",auxc1[3]+1),text[3],
                    rep(" ",auxc1[4]+1),text[4],"\n",sep="")
    auxc <- auxc+auxc1
    cat("  lambda: ",
        rep(ifelse(auxc[1] > aux1[1]," ",""),abs(auxc[1]-aux1[1])),sprintf("%.4f",x$lambda)," ",
        rep(ifelse(auxc[2] > aux2[1]," ",""),abs(auxc[2]-aux2[1])),sprintf("%.4f",x$lambda.std.error),
        "\n",sep="")
    cat("      xi: ",
        rep(ifelse(auxc[1] > aux1[2]," ",""),abs(auxc[1]-aux1[2])),sprintf("%.4f",x$xi)," ",
        rep(ifelse(auxc[2] > aux2[2]," ",""),abs(auxc[2]-aux2[2])),sprintf("%.4f",x$xi.std.error),
        "\n",sep="")
  
    if (length(x$beta) == 1) {
      i <- 1
      cat("    beta: ",
          rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
          rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.std.error[i])," ",
          rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),sprintf("%.4f",z.value[i])," ",
          rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),p.value2[i],
          ifelse(p.value[i] < 0.0001," ***",ifelse(p.value[i] < 0.01," **",
                 ifelse(p.value[i] < 0.05," *",ifelse(p.value[i] < 0.1," .","")))),
          "\n",sep="")
    } else {
      for (i in 1:length(x$beta)) {
        if (i-1 < 10) {
          cat("   beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.std.error[i])," ",
              rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),sprintf("%.4f",z.value[i])," ",
              rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),p.value2[i],
              ifelse(p.value[i] < 0.0001," ***",ifelse(p.value[i] < 0.01," **",
                     ifelse(p.value[i] < 0.05," *",ifelse(p.value[i] < 0.1," .","")))),
              "\n",sep="")
        } else {
          cat("  beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.std.error[i])," ",
              rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),sprintf("%.4f",z.value[i])," ",
              rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),p.value2[i],
              ifelse(p.value[i] < 0.0001," ***",ifelse(p.value[i] < 0.01," **",
                     ifelse(p.value[i] < 0.05," *",ifelse(p.value[i] < 0.1," .","")))),
              "\n",sep="")
        }
      }
    }
    cat("\n",sep="")
    cat("Summary statistic of the error for the TBS model\n",sep="")
    print(summary(x$error))
    cat("\n",sep="")
    if (length(x$beta) == 1) {
      cat("Estimated quantiles of time event\n",sep="")
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
  } else {
    cat("There was not obtained the convergence.\n",sep="")
  }
}

plot.tbs.survreg.mle <- function(x, ...) {
  if (x$convergence) {
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
      main <- c("Survival function (MLE)","Hazard function (MLE)","Error Distribution (MLE)")
    } else if (length(main) != 3) {
        stop("The 'main' length must be equal to 3")
    }
    axis.t <- seq(LB,UB,(UB-LB)/(h-1))

    par(ask=TRUE)
    for (k in 1:2) {
      if (attr(x$x,"plot") == 1) {
        if (k == 1) { ### Survival plot
          axis.y <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta,dist=x$error.dist)
          plot(axis.t,axis.y,type="l",xlim=c(0,UB),ylim=c(0,1),ylab=ylab[k],xlab=xlab[k],main=main[k], ...)
        } else { ### Hazard plot
          axis.y <- htbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta,dist=x$error.dist)
          plot(axis.t,axis.y,type="l",xlim=c(0,UB),ylab=ylab[k],xlab=xlab[k],main=main[k], ...)
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
        axis.y <- matrix(NA,length(axis.t),length(x$x))
        for (i in 1:length(x$x)) {
          if (attr(x$x,"plot") == 2) {
            if (k == 1) { ### Survival plot
              axis.y[,i] <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta*x$x[i],dist=x$error.dist)
            } else { ### Hazard plot
              axis.y[,i] <- htbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta*x$x[i],dist=x$error.dist)
            }
          } else {
            if (k == 1) { ### Survival plot
              axis.y[,i] <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,
                                   beta=(x$beta[1]+x$beta[2]*x$x[i]),dist=x$error.dist)
            } else { ### Hazard plot
              axis.y[,i] <- htbs(axis.t,lambda=x$lambda,xi=x$xi,
                                 beta=(x$beta[1]+x$beta[2]*x$x[i]),dist=x$error.dist)
            }
          }
          if (i == 1) {
            if (k == 1) { ### Survival plot
              plot(axis.t,axis.y[,i],xlim=c(0,UB),ylim=c(0,1),ylab=ylab[k],xlab=xlab[k],main=main[k],
                   type="l",lty=lty[i],lwd=lwd[i],col=col[i], ...)
            } else {
              plot(axis.t,axis.y[,i],xlim=c(0,UB),ylab=ylab[k],xlab=xlab[k],main=main[k],
                   type="l",lty=lty[i],lwd=lwd[i],col=col[i], ...)
            }
          } else {
            lines(axis.t,axis.y[,i],lty=lty[i],lwd=lwd[i],col=col[i])
          }
        }
      }
    }
  
    ## Histogram for error distribution:
    bound <- max(abs(c(min(x$error),max(x$error))))
    hist(x$error,probability=TRUE,xlim=c(-bound,bound),main=main[3],xlab=xlab[3],ylab=ylab[3])
    lines(seq(-bound,bound,2*bound/(h-1)),
          .choice(seq(-bound,bound,2*bound/(h-1)),xi=x$xi,dist=x$error.dist,type="d"))
  
    ## Q-Q plot for error distribution:
    qqplot(.choice(ppoints(length(x$error)), xi=x$xi, dist=x$error.dist, type="q"),x$error,
           main = expression("Q-Q plot for error"),xlab="Theoretical Quantiles",ylab="Sample Quantiles")
    qqline(x$error, distribution = function(p) .choice(p, xi=x$xi, dist=x$error.dist, type="q"),
           prob = c(0.25, 0.75), col = 2, lwd=2)
  
    par(ask=FALSE)
  } else {
    cat("There was not obtained the convergence.\n",sep="")
  }
}



