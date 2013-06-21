#!/usr/bin/env python

#ReactionMetParsingModule
import cPickle
import re


def RMParse(ModelRxnsFile,MapFlag): 
# Parse input files and return a dictionary mapping metabolites to reactions

	f1 = open(ModelRxnsFile,'r')
	Br = f1.read()
	f1.close()
	BrH = Br.split('\n')
	BesteReactants = []
	for item in BrH[:-1]:
		BesteReactants.append(item.split('\t'))
	
	BRrxn = []
	BRmet= []

	for item in BesteReactants:
		BRrxn.append(item[0])
		BRmet.append(item[1])

	# BesteRxns is the sorted, unique list of the reactions
	uBR = sorted(set(BRrxn))
	BesteMetsR = []
	for item in range(0,len(uBR)):
		BesteMetsR.append([])
		for rl,ml in zip(BRrxn,BRmet):
			if rl == uBR[item]:
				BesteMetsR[item].append(ml)

	BRmetsort = BesteMetsR
	for item in BRmetsort:
		item.sort()

	BRmetMap = BRmetsort
	#raise('WTF');
	if MapFlag  == 1:

		pklfile = open('BJdict.pkl','rb')
		dict1 = cPickle.load(pklfile)
		pklfile.close()
		for i1,i2 in zip(BRmetsort,BRmetMap):
			for m1 in range(0,len(i1)):
				if i1[m1] in dict1.keys():
					i2[m1] = dict1[i1[m1]]

	BRmetstr = []
	for item in BRmetMap:
		BRmetstr.append('-'.join(item))

	result = {}
	for met,rxn in zip(BRmetstr,uBR):
		try:
			result[met].append(rxn)
		except:
			result[met] = [rxn]

	return result

def  CommonMetOutputParse(CommonMetsList,MetDict1,MetDict2,CompareString): 
	#Converts the common metabolite list for 2 models into an informative output file

	## Outputing the results
	CommonMetsList = list(CommonMetsList)

	BRc = []
	FRc = []
	for item in CommonMetsList:
		BRc.append(', '.join(MetDict1[item]))
		FRc.append(', '.join(MetDict2[item]))

	filestr = CompareString[0] + CompareString[1] + 'commonMetrxns' + CompareString[2] + '.txt'
	
	fout = open(filestr,'w')
	fout.write('%s\t%s\t%s\n' %('Met',CompareString[0],CompareString[1]))

	for met,rB,rF in zip(CommonMetsList,BRc,FRc):
		fout.write('%s\t%s\t%s\n' %(met,rB,rF))
	fout.close()

def  CommonMetOutputParse3(CommonMetsList,MetDict1,MetDict2,MetDict3,CompareString): 
	#Converts the common metabolite list for 3 models into an informative output file

	## Outputing the results
	CommonMetsList = list(CommonMetsList)
	
	BcA = []
	FcA = []
	JcA = []
	for item in CommonMetsList:
		BcA.append(', '.join(MetDict1[item]))
		FcA.append(', '.join(MetDict2[item]))
		JcA.append(', '.join(MetDict3[item]))

	filestr = CompareString[0] + CompareString[1] + CompareString[2] + 'commonMetrxns' + CompareString[3] + '.txt'
	fout = open(filestr,'w')
	fout.write('%s\t%s\t%s\t%s\n' %('Met',CompareString[0],CompareString[1],CompareString[2]))
	for met,rB,rF,rJ in zip(CommonMetsList,BcA,FcA,JcA):
		fout.write('%s\t%s\t%s\t%s\n' %(met,rB,rF,rJ))
	fout.close()