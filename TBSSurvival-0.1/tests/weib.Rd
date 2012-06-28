% File TBSSurvival/tests/weib.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{weib}

%%%
\alias{dweib}
\alias{pweib}
\alias{qweib}
\alias{rweib}
\alias{hweib}
\alias{lweib}
\alias{bweib}

\title{The Weibull Distribution}
\description{
  Density, distribution function, quantile function, random
  generation, and hazard function for the Weibull distribution.
  The `lweib' is the log-likelihood function for survival given by (time,delta).
  The `bweib' performs the Bayesian estimation of the Weibull model with two parameters.
  The prior used is uniforme over real line for both parameters and the estimations are
  done by metroplois-hasting.
}
\usage{
 dweib(x,alpha,beta,mu=0)
 pweib(x,alpha,beta,mu=0)
 qweib(p,alpha,beta,mu=0)
 rweib(n,alpha,beta,mu=0)
 hweib(x,alpha,beta,mu=0)
 lweib(alpha,beta,mu=0,time,delta)
 bweib(guess,time,delta,burn,jump,size,scale)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If `length(n) > 1', the length is taken to be the number required.}
  \item{alpha}{scale parameter.}
  \item{beta}{shape parameter.}
  \item{mu}{location parameter. `mu = 0' is the Weibull of two parameters.}
  \item{time}{Survival/failure `time' array.}
  \item{delta}{Censoring indicator about `time' (0 means censored, 1 not censored).}
  \item{guess}{initial value of the Markov Chain.}
  \item{burn}{burn-in, number of firsts samples of posterior to not use.}
  \item{jump}{number of jump between each sample of posterior to avoid the problem of auto-correlation
              between the samples.}
  \item{size}{size of final sample of posterior.}
  \item{scale}{parameter of \code{metrop} function. Controls the acceptance rate.}
}
\value{
     `dweib' gives the density, `pweib' gives the distribution
     function, `qweib' gives the quantile function, and `rweib'
     generates random deviates. `hweib' is the hazard function.
     `lweib' is the log-likelihood function for the survival Weibull.
     `bweib' is the Bayesian estimation.
}
\examples{
  dweib(2,alpha=2,beta=1,mu=1)
}
