# Part of library 'TBSSurvival' for R (http://www.R-project.org)
#
# Copyright (C) 2010-2012, Adriano Polpo and D. Sinha.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
# 
# You should have received a copy of the GNU Library General Public
# License along with this library; A copy of the GNU General Public
# License is available on the World Wide Web at
# http://www.gnu.org/copyleft/gpl.html. You can also obtain it by
# writing to the Free Software Foundation, Inc., 51 Franklin St,
# Fifth Floor, Boston, MA  02110-1301, USA.

lik.weib <- function(par,time,delta)
{
  shape <- par[1]
  scale <- par[2]
  mu    <- 0
  if (length(par) == 3)
    mu <- par[3]

  out <- log(0)
  if ((scale > 0) && (min(time) > mu) && (shape > 0))
  {
    d.aux <- shape*(scale^(-shape))*((time[delta == 1]-mu)^(shape-1))*exp(-(((time[delta == 1]-mu)/scale)^shape))
    out <- sum(log(d.aux))    
    if (length(time[delta == 0]) != 0) {
      s.aux <- exp(-(((time[delta == 0]-mu)/scale)^shape))
      out <- out + sum(log(s.aux))
    }
  }
  return(out)
}

