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
tbs.survreg.mle <- function(formula,dist="norm",method=c("BFGS","Rsolnp","Nelder-Mead","CG","SANN"),verbose=FALSE,nstart=10,max.time=-1,ncore=1) {
  ## this meta-method only records the elapsed time and call the max likelihood estimation function
  ## for each of the methods given until one of them converges. It is supposed that at least one method
  ## is given, and that dist is one of those implemented by tbs.survreg.
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
    if(verbose) cat('No method has found a solution\n')
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

print.tbs.survreg.mle <- function(x, ...) {
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
  cat("-----------------------------------------------------\n",sep="")
  cat("TBS model with ",text.dist," error distribution.\n",sep="")
  cat("\n",sep="")
  cat("          point estimates (sd)\n",sep="")
  cat("  lambda: ",sprintf("%.4f",x$lambda),"  (",sprintf("%.4f",x$lambda.std.error),")\n",sep="")
  cat("      xi: ",sprintf("%.4f",x$xi),"  (",sprintf("%.4f",x$xi.std.error),")\n",sep="")
  for (i in 1:length(x$beta)) {
    cat("   beta",i-1,": ",sprintf("%.4f",x$beta[i]),"  (",sprintf("%.4f",x$beta.std.error[i]),")\n",sep="")
  }
  cat("\n",sep="")
  cat("    AIC: ",sprintf("%.4f",x$AIC),"\n",sep="")
  cat("   AICc: ",sprintf("%.4f",x$AICc),"\n",sep="")
  cat("    BIC: ",sprintf("%.4f",x$BIC),"\n",sep="")
  cat("\n",sep="")
  cat("Run time: ", sprintf("%.2f",x$run.time)," min \n",sep="")
  cat("-----------------------------------------------------\n",sep="")
  cat("\n",sep="")
}

summary.tbs.survreg.mle <- function(x, ...) {
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
  cat("-----------------------------------------------------\n",sep="")
  cat("TBS model with ",text.dist," error distribution.\n",sep="")
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
  for (i in 1:length(x$beta)) {
    cat("   beta",i-1,": ",
        rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
        rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.std.error[i])," ",
        rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),sprintf("%.4f",z.value[i])," ",
        rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),p.value2[i],
        ifelse(p.value[i] < 0.0001," ***",ifelse(p.value[i] < 0.01," **",
               ifelse(p.value[i] < 0.05," *",ifelse(p.value[i] < 0.1," .","")))),
        "\n",sep="")
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
  cat("-----------------------------------------------------\n",sep="")
  cat("\n",sep="")
}

plot.tbs.survreg.mle <- function(x, ...) {
  h <- 1000
  if (attr(x$x,"plot") == 1) {
    axis.t <- seq(0.0001,10,(10-0.0001)/(h-1))
    
  } else if (attr(x$x,"plot") == 2) {

  } else if (attr(x$x,"plot") == 3) {

  }
}


