\name{dweib}

%%%
\alias{dweib}
\alias{pweib}
\alias{qweib}
\alias{rweib}

\title{The Weibull Distribution}
\description{
  Density, distribution function, quantile function, and random
  generation for the Weibull distribution.
}
\usage{
 dweib(x,alpha,beta,mu=0)
 pweib(x,alpha,beta,mu=0)
 qweib(p,alpha,beta,mu=0)
 rweib(n,alpha,beta,mu=0)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If `length(n) > 1', the length is taken to be the number required.}
  \item{alpha}{scale parameter.}
  \item{beta}{shape parameter.}
  \item{mu}{location parameter. `mu = 0' is the Weibull of two parameters}
}
\value{
     `dweib' gives the density, `pweib' gives the distribution
     function, `qweib' gives the quantile function, and `rweib'
     generates random deviates.
}
\examples{
  dweib(2,alpha=2,beta=1,mu=1)
}
