#!/usr/bin/env python

import numpy as np

def quantilenorm(matrix):
	matrix = np.array(matrix)
	rankedmat = np.empty(shape=matrix.shape,dtype=int)
	for col in range(matrix.shape[1]):
		rankedmat[:,col] = ranking(matrix[:,col])

	sortedmat = np.empty(shape=matrix.shape)
	for col in range(matrix.shape[1]):
		sortedmat[:,col] = matrix[matrix[:,col].argsort(),col]

	rowavg = sortedmat.mean(axis=1)
	qnmat = np.empty(shape=matrix.shape,dtype=float)
	for count,avg in enumerate(rowavg):
		ix = np.where(rankedmat == count)
		ix = tuple(zip(*ix))
		for coord in ix:
			qnmat[coord[0],coord[1]] = rowavg[count]

	return qnmat

def ranking(vector):
	
	tmp = vector.argsort()
	ranks = np.empty(len(vector),int)
	ranks[tmp] = np.arange(len(vector))

	return ranks