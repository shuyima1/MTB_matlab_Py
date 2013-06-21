#!/usr/bin/env python

import re

##Beste Processing
file1 = open('BesteRxn.txt','r')
bw = file1.read()
file1.close()
BesteRxns = bw.split('\n')

file1 = open('BesteGPR.txt','r')
bw = file1.read()
file1.close()
BesteGPR = bw.split('\n')

BesteGenes = []
for item in BesteGPR:
	m = re.findall('Rv[0-9]{4}c?',item)
	BesteGenes.append(m)
	
	
# BesteDict = dict(zip(BesteRxns,BesteGenes))

BesteGenSort = BesteGenes
for item in BesteGenSort:
	item.sort()
	
BesteGenString = []
for item in BesteGenSort:
	BesteGenString.append('-'.join(item))
	
BesteGPRdict = {}

for gen,rxn in zip(BesteGenString,BesteRxns):
	try:
		BesteGPRdict[gen].append(rxn)
	except:
		BesteGPRdict[gen] = [rxn]
		
##Jamshidi Processing	
file1 = open('JamshidiRxn.txt','r')
bw = file1.read()
file1.close()
JamshidiRxns = bw.split('\n')

file1 = open('JamshidiGPR.txt','r')
bw = file1.read()
file1.close()
JamshidiGPR = bw.split('\n')

JamshidiGenes = []
for item in JamshidiGPR:
	m = re.findall('Rv[0-9]{4}c?',item)
	JamshidiGenes.append(m)
	
JamshidiGenSort = JamshidiGenes
for item in JamshidiGenSort:
	item.sort()
	
JamshidiGenString = []
for item in JamshidiGenSort:
	JamshidiGenString.append('-'.join(item))
	
JamshidiGPRdict = {}

for gen,rxn in zip(JamshidiGenString,JamshidiRxns):
	try:
		JamshidiGPRdict[gen].append(rxn)
	except:
		JamshidiGPRdict[gen] = [rxn]
		
##Fang Processing
file1 = open('FangRxn.txt','r')
bw = file1.read()
file1.close()
FangRxns = bw.split('\n')

file1 = open('FangGPR.txt','r')
bw = file1.read()
file1.close()
FangGPR = bw.split('\n')

FangGenes = []
for item in FangGPR:
	m = re.findall('Rv[0-9]{4}c?',item)
	FangGenes.append(m)

FangGenSort = FangGenes
for item in FangGenSort:
	item.sort()
	
FangGenString = []
for item in FangGenSort:
	FangGenString.append('-'.join(item))
	
FangGPRdict = {}

for gen,rxn in zip(FangGenString,FangRxns):
	try:
		FangGPRdict[gen].append(rxn)
	except:
		FangGPRdict[gen] = [rxn]
		
## Kalapanulak Processing
fK = open('Kalapanulak2009model.txt','r')
kw = fK.read()
fK.close()
fH = kw.split('\r')
KalInfo = []
for item in fH:
	KalInfo.append(item.split('\t'))
	
Kalrlong = []
Kalglong = []

for item in KalInfo[1:]:
	Kalrlong.append(item[0])
	Kalglong.append(item[4])

KalRxns = sorted(set(Kalrlong))
KalGenes = []
for item in range(0,len(KalRxns)):
	KalGenes.append([])
	for rl,gl in zip(Kalrlong,Kalglong):
		if rl == KalRxns[item]:
			KalGenes[item].append(gl)

KalGenSort = KalGenes
for item in KalGenSort:
	item.sort()

KalGenString = []
for item in KalGenSort:
	KalGenString.append('-'.join(item))

KalGPRdict = {}
for gen,rxn in zip(KalGenString,KalRxns):
	try:
		KalGPRdict[gen].append(rxn)
	except:
		KalGPRdict[gen] = [rxn]
			
			
## Comparing the Models
BesteKeys = set(BesteGPRdict.keys())
JamshidiKeys = set(JamshidiGPRdict.keys())
FangKeys = set(FangGPRdict.keys())

FJcommon = FangKeys.intersection(JamshidiKeys)
BFcommon = BesteKeys.intersection(FangKeys)
BJcommon = BesteKeys.intersection(JamshidiKeys)
BJFcommon = set.intersection(BesteKeys,JamshidiKeys,FangKeys)
BJFKcommon = set.intersection(BesteKeys,JamshidiKeys,FangKeys,KalKeys)
BKcommon= set.intersection(BesteKeys,KalKeys)
FKcommon = set.intersection(FangKeys,KalKeys)
JKcommon = set.intersection(JamshidiKeys,KalKeys)

## Outputing the results
BJFcgen = list(BJFcommon)

BRc = []
FRc = []
JRc = []
for item in BJFcommon:
	BRc.append(', '.join(BesteGPRdict[item]))
	FRc.append(', '.join(FangGPRdict[item]))
	JRc.append(', '.join(JamshidiGPRdict[item]))

fout = open('BJFcommonGPRrxns.txt','w')
fout.write('Genes\tBeste\tFang\tJamshidi\n')

for gen,rB,rF,rJ in zip(BJFcgen,BRc,FRc,JRc):
	fout.write('%s\t%s\t%s\t%s\n' %(gen,rB,rF,rJ))
fout.close()


BcA = []
FcA = []
JcA = []
KcA = []
for item in BJFKcommon:
	BcA.append(', '.join(BesteGPRdict[item]))
	FcA.append(', '.join(FangGPRdict[item]))
	JcA.append(', '.join(JamshidiGPRdict[item]))
	KcA.append(', '.join(KalGPRdict[item]))

fout = open('BJFKcommonGPRrxns.txt','w')
fout.write('Genes\tBeste\tFang\tJamshidi\tKalapanulak\n')
for gen,rB,rF,rJ,rK in zip(list(BJFKcommon),BcA,FcA,JcA,KcA):
	fout.write('%s\t%s\t%s\t%s\t%s\n' %(gen,rB,rF,rJ,rK)
fout.close()