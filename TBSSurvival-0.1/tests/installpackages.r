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
loc="/project/s221/cassio/Rlibs"
options(repos="http://stat.ethz.ch/CRAN/")
install.packages("mcmc",lib=loc)
install.packages("Rsolnp",lib=loc)
install.packages("normalp",lib=loc)
install.packages("eha",lib=loc)
install.packages("e1071",lib=loc)
install.packages("coda",lib=loc)
install.packages("truncnorm",lib=loc)
install.packages("R.utils",lib=loc)
