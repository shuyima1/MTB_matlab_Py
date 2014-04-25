#!/usr/bin/env python

import fileinput
import sys
from sklearn import preprocessing
from scipy.stats import ks_2samp
import numpy as np

data = []

count = -1
for line in fileinput.input():
	if fileinput.isfirstline():
		count = count + 1
		data.append([])
	data[count].append(line[:-1])

#print data

TRN = []
Exp = []
GenIDs = data[2]


for item in data[0]:
	TRN.append(item.split('\t'))

count = 0
for item in data[1]:
	Exp.append([])
	tmp = item.split(',')
	for element in tmp:
		Exp[count].append(float(element))
	count = count + 1

#print TRN[0]
#print Exp[0]
#print GenIDs[0]

del data

binarizer = preprocessing.Binarizer()
ExpBin = binarizer.transform(Exp)

for regpair in TRN:
	regix = GenIDs.index(regpair[0])
	tarix = GenIDs.index(regpair[1])

	tfonix = [i for i,x in enumerate(ExpBin[regix]) if x == 1]
	tfoffix = list(set(range(len(ExpBin[regix]))) - set(tfonix))