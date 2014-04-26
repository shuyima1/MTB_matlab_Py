#!/usr/bin/env python

from __future__ import division
import fileinput
import sys
from sklearn import preprocessing
from scipy.stats import ks_2samp
import numpy as np
import quantile


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

#count = 0
for count,item in enumerate(data[1]):
	Exp.append([])
	tmp = item.split(',')
	for element in tmp:
		Exp[count].append(float(element))
	#count = count + 1

#print TRN[0]
#print Exp[0]
#print GenIDs[0]

del data
Exp = np.array(Exp)
#binarizer = preprocessing.Binarizer()
#ExpBin = binarizer.transform(Exp)

#Qexp = quantile.quantilenorm(Exp)
ExpBin = []
#for count,row in enumerate(Qexp):
for count,row in enumerate(Exp):
	ExpBin.append([])
	for item in row:
		if item < 0:
			ExpBin[count].append(0)
		else:
			ExpBin[count].append(1)

tmp = []
for count,item in enumerate(ExpBin):
	tmp.append([])
	blah = list(item)
	for element in blah:
		tmp[count].append(str(int(element)))
print len(tmp)
print len(tmp[0])

tmp2 = []
for item in tmp:
	tmp2.append('\t'.join(item))

fout=open('PROMexpbin.txt','w')
for item in tmp2:
	fout.write('%s\n' %(item))

ExpBin = np.array(ExpBin)

probtfgene = [1.0]*len(TRN)
numerator = [1]*len(TRN)
denominator = [1]*len(TRN)
for count,regpair in enumerate(TRN):
	regix = GenIDs.index(regpair[0])
	tarix = GenIDs.index(regpair[1])

	tfonix = [i for i,x in enumerate(ExpBin[regix]) if x == 1]
	tfoffix = list(set(range(len(ExpBin[regix]))) - set(tfonix))

	if len(tfonix) > 0 and len(tfoffix) > 0:
		#ks = ks_2samp(Qexp[tarix][tfonix],Qexp[tarix][tfoffix])
		ks = ks_2samp(Exp[tarix][tfonix],Exp[tarix][tfoffix])

		if ks[1] < 0.05:
			numerator[count] = sum(ExpBin[tarix][tfoffix])
			denominator[count] = len(tfoffix)
			#probtfgene[count] = sum(ExpBin[tarix][tfoffix])/len(tfoffix)
			probtfgene[count] = numerator[count]/denominator[count]

print probtfgene[0]

fout = open('PROMprobtfgene.txt','w')
for p,n,d in zip(probtfgene,numerator,denominator):
	fout.write('%f\t%f\t%f\n' %(p,n,d))
fout.close()