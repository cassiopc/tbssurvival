% File TBSSurvival/man/dweib.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{tbs.survreg.mle}

%%%
\alias{tbs.survreg.mle}

\title{MLE of the TBS Model for Failure Data}
\description{
  This function perform the Maximum Likelihhod Estimation of the TBS model.
  The optimization is done by the function `optim' (or the package Rsolnp when available).
}
\usage{
 tbs.survreg.mle(formula,dist="norm",method=c("BFGS","Rsolnp","Nelder-Mead","CG","SANN"),verbose=FALSE,nstart=10,max.time=1,ncore=1)
}
\arguments{
  \item{formula}{A formula specification containing a Surv model with right-censored data as in the package survival.}
  \item{dist}{error distribution; \code{dist} = \code{"norm"}, \code{"doubexp"}, \code{"t"}, \code{"cauchy"} or \code{"logistic"}.}
  \item{method}{a vector of numerical methods to be used in the optimization. The function try all methods until reach the
                maximum with one of them or return `convergence = FALSE'.}
  \item{verbose}{Boolean to indicate the amount of output during the execution of the optimization.}
  \item{nstart}{Number of feasible initial points to guess when performing the optimization.}
  \item{max.time}{Maximum time (in minutes) to run the optimization.}
  \item{ncore}{If parallel computing is available, ncore tells the number of cores to use.}
}
\details{
  This function was developed to the problem without co-variables. Some changes maybe are need to work
  with co-variables (or maybe not).
}
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{std.error}{Standar error for MLE. If it is not possible to evaluate return \code{NA}}
  \item{log.lik}{The log-likelihood at parameters \code{par}.}
  \item{error.dist}{The error distribution chosen.}
  \item{AIC}{Akaike Information Criterion.}
  \item{AICc}{AICc is AIC with a second order correction for small sample sizes.}
  \item{BIC}{Bayesian Information Criterion.}
  \item{KS}{Komogorov-Simirnov Statistic. It is not evaluated with matrix \code{x} of co-variables.}
  \item{method}{Numerical method used to achive the MLE.}
  \item{convergence}{If \code{convergence} is \code{FALSE} then it was not possible to find the MLE.}
  \item{run.time}{Time spent with the function.}
}
\references{
  Meeker, W. and Escobar, L. (1998) \emph{Statistical Methods for Reliability Data}. Willey, ISBN 0-471-14328-6.
}
\seealso{
  \code{\link{tbs.survreg.bayes},\link{dtbs},\link{ptbs},\link{qtbs},\link{rtbs}}.
}
\examples{
# Alloy - T7987: data extracted from Meeker and Escobar (1998), pp. 131.
data(alloyT7987)
alloyT7987$time <- as.double(alloyT7987$time)
alloyT7987$delta <- as.double(alloyT7987$delta)

# MLE estimation with logistic error
formula <- Surv(alloyT7987$time,alloyT7987$delta == 1) ~ 1
tbs.mle <- tbs.survreg.mle(formula,dist="logistic")

# Kaplan-Meier estimation
km <- survfit(formula)

plot(km,ylab="",xlab="",xlim=c(min(alloyT7987$time),max(alloyT7987$time)),
     conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(alloyT7987$time),max(alloyT7987$time),
         (max(alloyT7987$time)-min(alloyT7987$time)-0.01)/1000)
title(ylab="R(t)",xlab="t: number of cycles (in thousands)",
      main="Reliability function (MLE)",
      sub="Alloy - T7987 (cf. Meeker and Escobar (1998), pp. 131)",
      cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=1,lwd.ticks=1,pos=0)
axis(2,lwd=1,lwd.ticks=1,pos=min(alloyT7987$time))
legend(200,0.95,c("Kaplan-Meier",
                  expression(textstyle(paste("TBS / ",sep="")) ~ epsilon
                             ~ textstyle(paste("~",sep="")) ~ Logistic)),
       col=c(1,2),lty=c(1,1),cex=1.1,lwd=c(1,2),bg="white")
lines(t,1-ptbs(t,lambda=tbs.mle$par[1],xi=tbs.mle$par[2],
               beta=tbs.mle$par[3],dist=tbs.mle$error.dist),type="l",
               lwd=2,col=2,lty=1)
}

