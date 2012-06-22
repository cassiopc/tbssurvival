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

tbs.survreg.mle <- function(formula,dist="norm",method=c("Rsolnp","BFGS","Nelder-Mead","CG","SANN"),verbose=FALSE) {
  ## this meta-method only records the elapsed time and call the max likelihood estimation function
  ## for each of the methods given until one of them converges. It is supposed that at least one method
  ## is given, and that dist is one of those implemented by tbs.survreg.
  initial.time <- .gettime()
  for(i in 1:length(method)) {
    out <- .tbs.survreg(formula,dist=dist,method=method[i],verbose=verbose)
    if(out$convergence)
      break
  }
  if (!out$convergence) {
    if(verbose) cat('No method achieved convergence\n')
    out$method <- NULL
  }
  ## run.time is returned in minutes
  out$call <- match.call()
  out$formula <- formula
  out$run.time <- .gettime() - initial.time
  return(out)
}
