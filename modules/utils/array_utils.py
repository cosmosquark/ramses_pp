import numpy as np

#Find the index of the closest value in an array
def argmin(array, value):
	idx = (np.abs(array-value)).argmin()
	return idx