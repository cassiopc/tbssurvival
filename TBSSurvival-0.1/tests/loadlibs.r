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

## This code is used for testing purposes. The TBSSurvival library does not
## depend on it for any of its functionalities

loadlibs <- function(libdir=NULL) {
  w <- options("warn")
  options("warn" = -1)
  if(require("TBSSurvival",quietly=TRUE)==FALSE) {
    library("survival")
    library("R.methodsS3",lib.loc=libdir)
    library("R.oo",lib.loc=libdir)
    library("R.utils",lib.loc=libdir)
    library("truncnorm",lib.loc=libdir)
    library("mcmc",lib.loc=libdir)
    library("normalp",lib.loc=libdir)
    library("eha",lib.loc=libdir)
    library("Rsolnp",lib.loc=libdir)
    library("e1071",lib.loc=libdir)
    library("coda",lib.loc=libdir)
    source("../R/tbs.survreg.bayes.r")
    source("../R/tbs.r")
    source("../R/tbs.survreg.mle.r")
    source("../R/local.r")
  }
  options("warn" = w[[1]])
}