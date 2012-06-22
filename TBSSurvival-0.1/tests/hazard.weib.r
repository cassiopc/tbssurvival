# Part of library 'GPTmodel' for R (http://www.R-project.org)
#
# Copyright (C) 2010, Adriano Polpo and D. Sinha.
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

hazard.weibull <- function(x,shape,scale)
  return((shape/scale)*((x/scale)^(shape-1)))

