# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2012-2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
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

## this function has the sole purpose of checking whether the arguments respect the
## needs of the other TBS functions' implementation. It also re-cast the arguments in
## case it is needed, but does not really perform calculations.
test.tbs <- function(lambda, xi, beta, x=NULL, time=NULL, type=NULL, p=NULL, n=NULL) {
  if (!is.numeric(xi))
    stop("xi is not a number")
  if (is.matrix(xi))
    stop("xi is matrix")
  if (xi <= 0)
    stop("xi <= 0")
  
  out   <- NULL
  out$x <- x
  out$beta <- beta

  if ((!is.numeric(lambda)) || (length(lambda) != 1))
    stop("lambda is not a number or length != 1")
  if (!is.numeric(beta))
    stop("beta is not a (vector) number")
  if (is.matrix(beta))
    stop("beta is matrix")
  if (!is.null(x)) {
    if (is.matrix(x)) {
      if (length(beta) != length(x[1,]))
        stop(paste("size of beta != ",length(x[1,]),sep=""))
    }
    else {
      if ((length(beta) != 1) && (length(beta) != length(x)))
        stop("size of beta is not conform")
    }
  }
  else {
    if (length(beta) > 1)
      stop("x is wrong or length(beta) > 1")  
  }
  if (lambda <= 0)
    stop("lambda <= 0")

  if (!is.null(type)) {
    if ((type == "d") || (type == "p")) {
      if (!is.numeric(time))
        stop("time is not a (vector) number")
      if (is.matrix(time))
        stop("time is matrix")
      if (any(time <= 0))
        stop("time <= 0")
      if (!is.null(x)) {
        if (is.matrix(x)) {
          if (length(time) != length(x[,1]))
            stop("length of time is different of length of x")
        }
        else {
          if (length(beta) == length(x)) {
            out$x <- matrix(x,1,length(x))
          } else {
            if (length(time) != length(x))
              stop("length of time is different of length of x")
            out$x <- matrix(x,length(x),1)
          }
        }
        out$beta <- matrix(beta,length(beta),1)
      }
      else {
        out$x <- matrix(1,length(time),1)
      }
    } else {
      if (type == "q") {
        if (!is.numeric(p))
          stop("p is not a (vector) number")
        if (is.matrix(p))
          stop("p is matrix")
        if (min(p) < 0)
          stop("p < 0")
        if (max(p) > 1)
          stop("p > 1")
      } else if (type == "r") {
          if (!is.numeric(n))
            stop("n is not a number")
          if (n %% 1 != 0)
            stop("n is not a integer number")
        }
      if (is.null(x)) {
        if (length(beta) > 1)
          stop("If x is omitted then beta must have length 1")
        out$x <- 1
      }
    }
  }

  return(out)
}
