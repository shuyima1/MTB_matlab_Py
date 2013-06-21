#!/usr/bin/env python

import re
import cPickle
import RxnMetParse

file1 = open('BJmap.txt','r')
bw = file1.read()
file1.close()
bH = bw.split('\n')
BJmap = []
for item in bH[1:448]:
	BJmap.append(item.split('\t'))

BJdict = dict(BJmap)

outfile = open('BJdict.pkl','wb')
cPickle.dump(BJdict,outfile)
outfile.close()

BRmetdict = RxnMetParse.RMparse('Beste_reactants.txt',1)
BRmetdict = RxnMetParse.RMparse('Beste_products.txt',1)
outfile = open('BesteMetdict.pkl','wb')
cPickle.dump(BRmetdict,outfile)
cPickle.dump(BPmetdict,outfile)
outfile.close()

JRmetdict = RxnMetParse.RMparse('Jamshidi_reactants.txt',0)
JPmetdict = RxnMetParse.RMparse('Jamshidi_products.txt',0)
outfile = open('JamshidiMetdict.pkl','wb')
cPickle.dump(JRmetdict,outfile)
cPickle.dump(JPmetdict,outfile)
outfile.close()

FRmetdict = RxnMetParse.RMparse('Fang_reactants.txt',0)
FPmetdict = RxnMetParse.RMparse('Fang_products.txt',0)
outfile = open('FangMetdict.pkl','wb')
cPickle.dump(FRmetdict,outfile)
cPickle.dump(FPmetdict,outfile)
outfile.close()

BRmetKeys = set(BRmetdict.keys())
BPmetKeys = set(BPmetdict.keys())
JRmetKeys = set(JRmetdict.keys())
JPmetKeys = set(JPmetdict.keys())
FRmetKeys = set(FRmetdict.keys())
FPmetKeys = set(FPmetdict.keys())

FJrcommon = FRmetKeys.intersection(JRmetKeys)
FJpcommon = FangKeys.intersection(JamshidiKeys)
RxnMetParse.CommonMetOutputParse(FJrcommon,FRmetdict,JRmetdict,'FJR')
RxnMetParse.CommonMetOutputParse(FJpcommon,FPmetdict,JPmetdict,'FJP')

BFrcommon = BRmetKeys.intersection(FRmetKeys)
BFpcommon = BesteKeys.intersection(FangKeys)
RxnMetParse.CommonMetOutputParse(BFrcommon,BRmetdict,FRmetdict,'BFR')
RxnMetParse.CommonMetOutputParse(BFpcommon,BPmetdict,FPmetdict,'BFP')

BJrcommon = BRmetKeys.intersection(JRmetKeys)
BJpcommon = BesteKeys.intersection(JamshidiKeys)
RxnMetParse.CommonMetOutputParse(BJrcommon,BRmetdict,JRmetdict,'BJR')
RxnMetParse.CommonMetOutputParse(BJpcommon,BPmetdict,JPmetdict,'BJP')

BJFrcommon = set.intersection(BRmetKeys,JRmetKeys,FRmetKeys)
BJFpcommon = set.intersection(BesteKeys,JamshidiKeys,FangKeys)
RxnMetParse.CommonMetOutputParse3(BJrcommon,BRmetdict,FRmetdict,JRmetdict,'BFJR')
RxnMetParse.CommonMetOutputParse3(BJpcommon,BPmetdict,FPmetdict,JPmetdict,'BFJP')

