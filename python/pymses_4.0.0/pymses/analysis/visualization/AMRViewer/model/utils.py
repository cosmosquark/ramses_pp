# License:
#   Copyright (C) 2011 Thomas GUILLET, Damien CHAPON, Marc LABADENS. All Rights Reserved.
#
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.
from numpy import array, log10, argmin, abs

def get_appropriate_unit(ph_value, unit_dict, nearest_log10=1.0):
	"""
	"""
	c = array([ph_value.express(unit_dict[k]) for k in unit_dict.keys()])
	if c.any():
		lv = abs(log10(c)-nearest_log10)
	else:
		lv =c
	ilmin = argmin(lv)
	unit = unit_dict.keys()[ilmin]
	ph_v = c[ilmin]
	return (ph_v, unit)
	
