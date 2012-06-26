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
## hazard function for the TBS
htbs <- function(time,lambda=1,xi=1,beta=1,x=NULL,dist="norm") {
  return(dtbs(time,lambda,xi,beta,x,dist)/(1-ptbs(time,lambda,xi,beta,x,dist)))
}

#######################################################################
## density function for the TBS
dtbs <- function(time,lambda=1,xi=1,beta=1,x=NULL,dist="norm") {
  aux <- .test.tbs(lambda,xi,beta,x,dist,time=time,type="d")
  out <- ((1/time)*(abs(log(time))^(lambda-1))*
          .choice(c(.g.lambda(log(time),lambda)-.g.lambda(c(aux$x%*%aux$beta),lambda)),xi,dist,"d"))
  return(c(out))
}

#######################################################################
## distribution function for the TBS
ptbs <- function(time,lambda=1,xi=1,beta=1,x=NULL,dist="norm") {
  aux <- .test.tbs(lambda,xi,beta,x,dist,time=time,type="d")
  out <- .choice(c(.g.lambda(log(time),lambda)-.g.lambda(c(aux$x%*%aux$beta),lambda)),xi,dist,"p")
  return(c(out))
}

#######################################################################
## quantile function for the TBS
qtbs <- function(p,lambda=1,xi=1,beta=1,x=NULL,dist="norm") {
  aux <- .test.tbs(lambda,xi,beta,x,dist,type="q",p=p)

  aux2 <- c(aux$x%*%aux$beta)
  if (length(aux2) == 1)
    out <- c(exp(.g.lambda.inv(.g.lambda(aux2,lambda)+
             .choice(p,xi,dist,"q"),lambda)))
  else {
    out <- matrix(0,length(p),length(aux2))
    for (i in 1:length(aux2))
      out[,i] <- c(exp(.g.lambda.inv(.g.lambda(aux2[i],lambda)+.choice(p,xi,dist,"q"),lambda)))
  }

  return(out)
}

#######################################################################
## random generation for the TBS
rtbs <- function(n,lambda=1,xi=1,beta=1,x=NULL,dist="norm") {
 aux <- .test.tbs(lambda,xi,beta,x,dist,type="r",n=n)

 aux2 <- c(aux$x%*%aux$beta)
 erro <- .choice(n,xi,dist,"r")
 out  <- exp(.g.lambda.inv(.g.lambda(aux2,lambda) + erro,lambda))
 ## just to avoid time zero:
 out  <- ifelse(out == 0,10e-100,out)

 return(out)
}
