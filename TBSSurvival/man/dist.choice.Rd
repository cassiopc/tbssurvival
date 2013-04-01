% File TBSSurvival/man/dist.choice.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2013, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{dist.choice}

%%%
\alias{test.tbs}

\title{The Wrapper for new functions to use with TBS}
\description{
  Density, distribution function, quantile function, random
  generation and hazard function for the Transform-Both-Sides (TBS) distribution.
}
\usage{
  dist.choice(dist="norm")
}
\arguments{
  \item{dist}{Distribution of error, dist = "norm", "t", "doubexp", "cauchy" or "logistic". A choice "all" can also be
given, in which case a list with all the available distributions will be generated.}
}
\details{
  This wrapper creates a list with the necessary information about a distribution to be used with Transform-Both-Sides (TBS) model.
  It returns a list with six items, namely the density function, distribution function, quantile function, random generation function, test
  function, and name (a string). For example, with dist = "norm" it gives list(dnorm,pnorm,qnorm,rnorm,test.tbs,"norm"). The idea is that
  the user can implement its own list to be used with the TBS, as long as it is a unimodal symmetric distribution, and replace the call
  of dist.choice with their own list. The \link{test.tbs} function is used to check the arguments of the other functions. You may use the
  TBS provided function for that.
}
\value{
     A list(d.dist,p.dist,q.dist,r.dist,test.tbs,name.dist) according with the chosen distribution.
     The standard function for testing is test.tbs(lambda, xi, beta, x=NULL, time=NULL, type=NULL, p=NULL, n=NULL)
     The others (d.dist,p.dist,q.dist,r.dist) have to accept exactly two arguments: the first is the actual argument to the function, and the second is a parameter of the it. 
}
\examples{
  ## this will return list(dcauchy,pcauchy,qcauchy,rcauchy,test.tbs,"cauchy")
  dist.choice("cauchy")
}
