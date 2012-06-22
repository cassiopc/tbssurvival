% File TBSSurvival/man/bayes.weib.est.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 2 or later

\name{bayes.weib.est}

%%%
\alias{bayes.weib.est}

\title{Bayesian Estimation of the Weibull Model for Failure Data}
\description{
  This function perform the Bayesian estimation of the Weibull model with two parameters.
  The prior used is uniforme over real line for both parameters and the estimations are
  done by metroplois-hasting (using the function \code{metrop} availible at package \code{mcmc}).
}
\usage{
 bayes.weib.est(kick,time,delta,burn,jump,size,scale)
}
\arguments{
  \item{kick}{initial value of the Markov Chain.}
  \item{time}{vector of observed failure times.}
  \item{delta}{vector that indicates if the time is a censor (delta = 0) or not (delta = 1)}
  \item{burn}{burn-in, number of firsts samples of posterior to not use.}
  \item{jump}{number of jump between each sample of posterior to avoid the problem of auto-correlation
              between the samples.}
  \item{size}{size of final sample of posterior.}
  \item{scale}{parameter of \code{metrop} function. Controls the acceptance rate.}
}
