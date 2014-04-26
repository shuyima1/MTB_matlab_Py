#!/usr/bin/env python

import numpy as np

def quantilenorm(matrix):
	matrix = np.array(matrix)
	rankedmat = np.empty(shape=matrix.shape,dtype=int)
	for col in range(matrix.shape[1]):
		rankedmat[:,col] = ranking(matrix[:,col])

	return qnmat

def ranking(vector):
	
	tmp = vector.argsort()
	ranks = np.empty(len(vector),int)
	ranks[tmp] = np.arange(len(vector))

	return ranks