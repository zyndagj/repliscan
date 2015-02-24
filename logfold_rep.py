#!/usr/bin/env python

import os,sys
import math
import argparse
import numpy as np
from multiprocessing import Pool, cpu_count

def run_logFC(fList):
	beds = map(lambda y: y[1]+".bedgraph", fList)
	print len(beds), beds
	vals = []
	print "Parsing Files"
	locations = processFiles(vals, beds)
	vals = np.array(vals)
	normVals = normalize(vals)
	del vals
	for bIndex in xrange(1,len(beds)):
		outName = fList[bIndex][1]+"_logFC.bedgraph"
		if not os.path.exists(outName):
			print "Making %s"%(outName)
			OF = open(outName,'w')
			subV = (normVals[bIndex,:]+1)/(normVals[0,:]+1)
			lVals = np.log2(subV)
			for i in xrange(len(locations)):
				outStr = '\t'.join(map(str, locations[i]+[lVals[i]]))
				OF.write(outStr+'\n')
			OF.close()

def smooth(level, fList, chromDict):
	beds = map(lambda y: y[1]+"_logFC.bedgraph", fList)
	print "Smoothing:",beds
	for bed in beds:
		outFile = '%s_%i.smooth.bedgraph'%(bed.rstrip('.bedgraph'),level)
		os.system('rm %s'%(outFile))
		for chr in chromDict:
			os.system('grep "^%s\t" %s | cut -f 4 > tmp.val' % (chr, bed))
			os.system('grep "^%s\t" %s | cut -f 1-3 > tmp.loc'%(chr, bed))
			os.system('./wavelets --level %i --to-stdout --boundary reflected --filter Haar tmp.val > tmp.smooth'%(level))
			os.system('paste tmp.loc tmp.smooth >> %s'%(outFile))
	os.system('rm tmp.val tmp.loc tmp.smooth')

def geometricMean(M):
	'''
	Returns the geometric mean of a numpy array.

	Parameters
	============================================
	M	numpy array (matrix)
	
	>>> geometricMean(np.array([[1,2,3],[4,5,6]]))
	array([ 2.		,  3.16227766,  4.24264069])
	>>> geometricMean(np.array([[0.1, 1, 2],[10,1,0.2]]))
	array([ 1.		,  1.		,  0.63245553])
	'''
	return np.prod(M,axis=0)**(1.0/M.shape[0])

def sizeFactors(npVals):
	'''
	Calculates size factors like DESeq. http://genomebiology.com/2010/11/10/R106
	- equation 5

	Parameters
	============================================
	M	numpy array (matrix)

	>>> np.round(sizeFactors(np.array([[0.1,1,2],[10,1,0.2]])),2)
	array([ 1.,  1.])
	>>> sizeFactors(np.array([[-1,1],[0,3]]))
	Traceback (most recent call last):
	 ...
	ArithmeticError: Negative values in matrix
	>>> np.round(sizeFactors(np.array([[1,2,3],[4,5,6]])),2)
	array([ 0.63,  1.58])
	'''
	if np.any(npVals < 0):
		raise ArithmeticError("Negative values in matrix")
	piArray = geometricMean(npVals)
	gZero = piArray > 0
	return np.median(npVals[:,gZero]/piArray[gZero], axis=1)

def normalize(vals):
	'''
	Normalize the counts using DESeq's method

	>>> np.round(normalize(np.array([[1,2,3],[4,5,6]])),2)
	array([[ 1.58,  3.16,  4.74],
		   [ 2.53,  3.16,  3.79]])
	'''
	sf = sizeFactors(vals)
	return np.array(vals/np.matrix(sf).T)

def processFiles(vals, files):
	p = Pool(cpu_count())
	for l,v in p.map(parseBG, files):
		vals.append(v)
	return l

def parseBG(inFile):
	lines = open(inFile,'r').readlines()
	tmp = map(lambda y: y.rstrip('\n').split('\t'), lines)
	locs = map(lambda y: y[:3], tmp)
	vals = np.array(map(lambda y: np.float(y[3]), tmp))
	vals[vals < 2] = 0
	return locs, vals

def readFAI(inFile):
	'''
	Returns (length, seek, read length, readlength + newline)
	'''
	#Pt     140384  50      60      61
	chromDict = {}
	for line in open(inFile,'r'):
		tmp = line.split('\t')
		chromDict[tmp[0]] = tuple(map(int, tmp[1:]))
	return chromDict

def parseIN(inFile):
	fList = []
	for line in open(inFile,'r'):
		tmp = line.rstrip('\n').split('\t')
		fList.append(tuple(tmp))
	return fList

def makeBedgraph(fList, fasta, size):
	fai = fasta+".fai"
	bed = fasta+".bed"
	if not os.path.exists(fai):
		print "Making fai index"
		os.system("samtools faidx %s"%(fasta))
	if not os.path.exists(bed):
		print "Making 500bp window bed"
		os.system("bedtools makewindows -g %s -w %i | sort -S 10G -k1,1 -k2,2n > %s"%(fai,size,bed))
	for bam, prefix in fList:
		bg = prefix+'.bedgraph'
		if not os.path.exists(bg):
			print "Generating intersect for %s"%(bam)
			os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%(bed,bam,bg))

def main():
	parser = argparse.ArgumentParser(description="Performs log-fold analysis on bam files. The main input is a text file that contains a list of bam files. The first line is used as the control, or input.")
	parser.add_argument("infile", metavar="FILE", help="File with list of bams")
	parser.add_argument("-F",metavar='FASTA',help="Fasta file", required=True)
	parser.add_argument("-L",metavar='INT', help="Smoothing level", default=2, type=int)
	parser.add_argument("-S",metavar='INT', help="Bin size (Default: %(default)s)", default=500, type=int)
	args = parser.parse_args()
	if os.path.splitext(args.F)[1] in ['.fasta','.fa']:
		fai = args.F+'.fai'
	else:
		sys.exit("Please specify a fasta file\n")
	fList = parseIN(args.infile)
	makeBedgraph(fList, args.F, args.S)
	chromDict = readFAI(fai)
	run_logFC(fList)
	sys.exit()
	for i in range(1,args.L+1):
		smooth(i, fList, chromDict)

if __name__ == "__main__":
	main()
