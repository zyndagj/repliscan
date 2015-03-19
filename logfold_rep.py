#!/usr/bin/env python

import os,sys
import math
import argparse
import numpy as np
from multiprocessing import Pool, cpu_count
from glob import glob
import gc
from copy import deepcopy
from array import array

myColors = ("#85BEFF", "#986300", "#009863", "#F2EC00", "#F23600", "#C21BFF", "#85FFC7")
colorDict = {frozenset([0]):myColors[0], frozenset([1]):myColors[1], frozenset([2]):myColors[2], frozenset([2,1]):myColors[3], frozenset([2,0]):myColors[4], frozenset([1,0]):myColors[5], frozenset([2,1,0]):myColors[6]}
zero = 0.0004

def printMem():
	mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_AVPHYS_PAGES')  # e.g. 4015976448
	mem_gib = mem_bytes/(1024.**3)  # e.g. 3.74
	print "Memory:",mem_gib

def run_logFC(fList, thresh=2.5):
	beds = map(lambda y: y[1]+".bedgraph", fList)
	if not os.path.exists(fList[1][1]+"_logFC.bedgraph"):
		print "Making LogFold Files From:", beds
		L, naVals = processFiles(beds)
		gc.collect()
		threshVals(naVals[0,:], thresh) # raise low counts to minimal coverage in G1
		normVals = normalize(naVals)
		for bIndex in xrange(1,len(beds)):
			outName = fList[bIndex][1]+"_logFC.bedgraph"
			if not os.path.exists(outName):
				print "Making %s"%(outName)
				OF = open(outName,'w')
				subV = (normVals[bIndex,:]+1.0)/(normVals[0,:]+1.0)
				lVals = np.log2(subV)
				for i in xrange(len(L[0])):
					outStr = '%s\t%i\t%i\t%.4f\n' % (L[0][i], L[1][i], L[2][i], lVals[i])
					#outStr = '\t'.join(map(str, locations[i]+[lVals[i]]))
					OF.write(outStr)
				OF.close()
		del subV, lVals, L, normVals
	gc.collect()

def threshVals(vals, thresh):
	lgz = vals > 0
	minVal = np.percentile(vals[lgz], thresh)
	lmin = np.logical_and(vals < minVal, lgz)
	vals[lmin] = minVal

def normalize(vals):
	'''
	Normalize the counts using DESeq's method

	>>> np.round(normalize(np.array([[1,2,3],[4,5,6]])),2)
	array([[ 1.58,  3.16,  4.74],
		   [ 2.53,  3.16,  3.79]])
	'''
	sf = sizeFactors(vals)
	nA = np.array(vals/np.matrix(sf).T)
	return nA

def processFiles(files):
	p = Pool(12)
	vals = []
	for i in p.imap(parseVals, files):
		vals.append(i)
	p.close()
	p.join()
	del p
	gc.collect()
	return (parseLocs(files[0]), np.array(vals,dtype=np.float32))

def parseLocs(inFile):
	chroms = []
	starts = array('I',[])
	ends = array('I',[])
	for line in open(inFile,'r'):
		tmp = line.split('\t')
		chroms.append(tmp[0])
		starts.append(int(tmp[1]))
		ends.append(int(tmp[2]))
	return (chroms, starts, ends)

def parseVals(inFile):
	vals = array('f',[])
	for line in open(inFile,'r'):
		tmp = line.rstrip('\n').split('\t')
		vals.append(float(tmp[3]))
	return vals

def smooth(level, fList, chromDict):
	sortedChroms = sorted(chromDict.keys()[:])
	beds = map(lambda y: y[1]+"_logFC.bedgraph", fList[1:])
	print "Smoothing:",beds
	for bed in beds:
		outFile = '%s_%i.smooth.bedgraph'%(bed.rstrip('.bedgraph'),level)
		if not os.path.exists(outFile):
			for chr in sortedChroms:
				os.system('grep "^%s\t" %s | cut -f 4 > tmp.val' % (chr, bed))
				os.system('grep "^%s\t" %s | cut -f 1-3 > tmp.loc'%(chr, bed))
				os.system('wavelets --level %i --to-stdout --boundary reflected --filter Haar tmp.val > tmp.smooth'%(level))
				os.system('paste tmp.loc tmp.smooth >> %s'%(outFile))
				os.system('rm tmp.val tmp.loc tmp.smooth')
	gc.collect()

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
	sf = np.median(npVals[:,gZero]/piArray[gZero], axis=1)
	return sf

def readFAI(inFile):
	'''
	Returns (length, seek, read length, readlength + newline)
	'''
	#Pt     140384  50      60      61
	chromDict = {}
	for line in open(inFile,'r'):
		tmp = line.split('\t')
		#chromDict[tmp[0]] = tuple(map(int, tmp[1:]))
		chromDict[tmp[0]] = int(tmp[1])
	return chromDict

def parseIN(inFile):
	fList = []
	for line in open(inFile,'r'):
		tmp = line.rstrip('\n').split('\t')
		fList.append(tuple(tmp))
	return fList

def makeBedgraph(fList, fasta, size):
	fai = fasta+".fai"
	bed = "%s.%i.bed"%(fasta,size)
	if not os.path.exists(fai):
		print "Making fai index"
		os.system("samtools faidx %s"%(fasta))
	else: print "%s exists already"%(fai)
	if not os.path.exists(bed):
		print "Making %ibp window bed"%(size)
		os.system("bedtools makewindows -g %s -w %i | sort -S 20G -k1,1 -k2,2n > %s"%(fai,size,bed))
	else: print "%s exists already"%(bed)
	for bam, prefix in fList:
		bg = prefix+'.bedgraph'
		if not os.path.exists(bg):
			print "Generating intersect for %s"%(bam)
			os.system("bedtools bamtobed -i %s | cut -f 1-3 | sort -S 20G -k1,1 -k2,2n | bedtools intersect -a %s -b stdin -c -sorted > %s"%(bam,bed,bg))

def main():
	parser = argparse.ArgumentParser(description="Performs log-fold analysis on bam files. The main input is a text file that contains a list of bam files. The first line is used as the control, or input.")
	parser.add_argument("infile", metavar="FILE", help="File with list of bams")
	parser.add_argument("-F",metavar='FASTA',help="Fasta file", required=True)
	parser.add_argument("-L",metavar='INT', help="Smoothing level (Default: %(default)s)", default=2, type=int)
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
	smooth(args.L, fList, chromDict)
	gc.collect()
	makeGFF(fList, chromDict, args.L, args.S)
	print "Done"

def parseLocations(chroms):
	ctmp = ""
	start = 0
	locDict = {}
	for i in xrange(len(chroms)):
		chrom = chroms[i]
		if chrom != ctmp:
			if ctmp != "":
				locDict[ctmp] = (start, i)
			ctmp = chrom
			start = i
	locDict[ctmp] = (start,i+1)
	return locDict

def makeGFF(fList, chromDict, level, S):
	sortedChroms = sorted(chromDict.keys()[:])
	beds = map(lambda y: "%s_logFC_%i.smooth.bedgraph"%(y[1],level), fList[1:])
	names = map(lambda y: y[1], fList[1:])
	for name in names:
		outGFF = '%s_logFC_%i.smooth.gff3'%(name,level)
		open(outGFF,'w').write('##gff-version 3\n#track name="%s LogFold %ibp" gffTags=on\n' % (name,S))
		open("logFC_segmentation.gff3",'w').write('##gff-version 3\n#track name="Segmentation %ibp" gffTags=on\n'%(S))
	print "Parsing :",beds
	L, vals = processFiles(beds)
	locDict = parseLocations(L[0])
	for chrom in sortedChroms:
		s,e = locDict[chrom]
		tmpChr = L[0][s:e]
		tmpS = L[1][s:e]
		tmpE = L[2][s:e]
		maskM = np.zeros((len(names),len(tmpChr)),dtype=np.bool)
		for i in xrange(len(beds)):
			outGFF = '%s_logFC_%i.smooth.gff3'%(names[i],level)
			outSignal = vals[i,s:e]
			bMask = np.zeros(len(outSignal), dtype=np.bool)
			bMask[outSignal > zero] = 1
			maskM[i,:]=bMask[:]
			rBounds = calcRegionBounds(bMask)
			OF = open(outGFF,'a')
			count = 1
			for rS,rE in rBounds:
				OF.write("%s\t.\tgene\t%i\t%i\t.\t.\t.\tID=gene%i;color=%s;\n" % (tmpChr[rS], tmpS[rS]+1, tmpE[rE], count, myColors[i]))
				count += 1
			OF.close()
		OS = open('logFC_segmentation.gff3','a')
		count = 1
		setSigI = set(range(len(names)))
		for s in powerSet(range(len(names))):
			setS = set(s)
			tA = np.ones(maskM.shape[1], dtype=np.bool)
			for i in s: #intersection
				tA = np.logical_and(maskM[i,:], tA)
			for i in list(setSigI-setS): #remove differences
				tA[maskM[i,:]] = 0
			rBounds = calcRegionBounds(tA)
			name = ''.join(map(lambda y: names[y],sorted(s)))
			for rS,rE in rBounds:
				OS.write("%s\t.\tgene\t%i\t%s\t.\t.\t.\tID=gene%i;Name=%s;color=%s;\n" % (tmpChr[rS], tmpS[rS]+1, tmpE[rE], count, name, colorDict[frozenset(s)]))
				count += 1
		OS.close()

def powerSet(L):
	'''
	Removes the empty list from powerset

	>>> powerSet([1,2])
	[[1], [0], [1, 0]]
	'''
	def powerHelp(A):
		'''
		Builds powerset recursively.

		>>> powerHelp([1,2])
		[[], [1], [0], [1, 0]]
		'''
		if not A:
			return [[]]
		ret = powerHelp(A[1:])
		return ret+[i+[A[0]] for i in ret]
	return powerHelp(L)[1:]

def haarVec(vals, p=80):
	hC = pywt.wavedec(vals,'haar')
	cutVal = np.percentile(np.abs(np.concatenate(hC)), p)
	for A in hC:
		A[np.abs(A) < cutVal] = 0
	tVals = pywt.waverec(hC,'haar')
	return tVals[:len(vals)]

def hannVec(vals, windowSize=20):
	# Smooths using a hanning window
	w = np.hanning(windowSize)
	return np.convolve(w/w.sum(), vals, mode='same')

def mergeRegions(counter, distThresh=0):
	'''
	Merges regions that are closer than (<) the distance threshold.

	Parameters
	=============================
	counter		Binary counter array
	distThresh	Max distance threshold

	>>> A=np.array([1,1,0,0,0,1,1])
	>>> mergeRegions(A)
	>>> A
	array([1, 1, 0, 0, 0, 1, 1])
	>>> mergeRegions(A, distThresh=3)
	>>> A
	array([1, 1, 0, 0, 0, 1, 1])
	>>> mergeRegions(A, distThresh=4)
	>>> A
	array([1, 1, 1, 1, 1, 1, 1])
	'''
	bounds = calcRegionBounds(counter)
	for i in xrange(len(bounds)-1):
		start0, end0 = bounds[i]
		start1, end1 = bounds[i+1]
		if start1-end0-1 < distThresh:
			counter[start0:end1] = 1
def calcRegionBounds(counter):
	'''
	Returns the new lower and upper bounds over overlapped regions.

	Parameters
	=============================
	counter		Binary counter array

	>>> calcRegionBounds(np.array([1,1,0,0,1,1,1,0,0,1,1]))
	array([[ 0,  1],
	       [ 4,  6],
	       [ 9, 10]])
	'''
	d = np.diff(counter)
	idx, = d.nonzero()
	if counter[0]:
		idx = np.r_[-1, idx]
	if counter[-1]:
		idx = np.r_[idx, counter.size-1]
	idx.shape = (-1,2)
	idx[:,0] += 1
	return idx

if __name__ == "__main__":
	main()
