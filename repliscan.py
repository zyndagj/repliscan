#!/usr/bin/env python

import os, sys, math, argparse, gc, colorsys
import numpy as np
from multiprocessing import Pool, cpu_count, Process, Pipe
import subprocess as sp
from glob import glob
from copy import deepcopy
from array import array
import scipy.optimize
import scipy.interpolate
import scipy.misc
import scipy.stats as ss
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

myColors = ("#2250F1","#1A8A12","#FB0018","#FFFD33","#EA3CF2","#28C5CC","#FAB427")
colorDict = {frozenset([0]):myColors[0], frozenset([1]):myColors[1], frozenset([2]):myColors[2], frozenset([2,1]):myColors[3], frozenset([2,0]):myColors[4], frozenset([1,0]):myColors[5], frozenset([2,1,0]):myColors[6]}

class argChecker():
	def __init__(self, options, afterValid):
		self.options = options
		self.av = afterValid
	def check(self, x):
		if x in self.options:
			return x
		else:
			raise argparse.ArgumentTypeError("%s not a valid %s"%(x, self.av))

def main():
	parser = argparse.ArgumentParser(description="Performs log-fold analysis on bam files.",formatter_class=argparse.RawDescriptionHelpFormatter, epilog='''\
Input TXT File:
  Each line of the text file needs to contain
  a short name describing the sample and then a
  list of bam files corresponding to that name,
  all separated by tabs.

  The first line of this
  file needs to be the control (G1).
  All subsequent lines need to be listed
  sequentially according to experimental time.

  Example TXT File:
  G1	G1_001.bam	G1_002.bam
  ES	ES_001.bam
  MS	MS_001.bam	MS_L1.bam	MS_L2.bam

Methods to handle replicates:
  - sum (Default)
  - median
  - mean
  - min
  - max''')
	parser.add_argument("infile", metavar="FILE", help="File with list of bams")
	parser.add_argument("-F",metavar='FASTA',help="Fasta file", required=True)
	parser.add_argument("-L",metavar='INT', help=\
		"Smoothing level (Default: %(default)s)", default=3, type=int)
	parser.add_argument("-S",metavar='INT', help=\
		"Bin size (Default: %(default)s)", default=1000, type=int)
	parser.add_argument("-C",metavar='STR', help=\
		"Replicate reduction method (Default: %(default)s)", default="sum",\
		type=argChecker(('sum','median','mean','min','max'), 'reducer').check)
	parser.add_argument("--use",metavar='STR', help=\
		"Data to use for smoothing/segmentation (log|ratio Default: %(default)s)",\
		default="ratio", type=argChecker(('log','ratio'),'format').check)
	parser.add_argument("--norm", metavar='STR', help=\
		"Normalization Method (DESeq|Coverage) (Default: %(default)s)", default="Coverage",\
		type=argChecker(('DESeq','Coverage'),'normalization method').check)
	parser.add_argument("--rep", metavar='STR', help=\
		"Replication Method (threshold|auto|percent) (Default: %(default)s)", default="auto",\
		type=argChecker(('threshold','auto','percent'),'replication method').check)
	parser.add_argument("-T", metavar='Float', help=\
		"Threshold Level (Default: %(default)s)", default=0.0, type=float)
	parser.add_argument("-P", metavar='Float', help=\
		"Percent Cut (Default: %(default)s)", default=2.0, type=float)
	parser.add_argument("--seg", metavar='STR', help=\
		"Segmentation Method (binary|proportion) (Default: %(default)s)", default="proportion",\
		type=argChecker(('binary','proportion'),'segmentation method').check)
	parser.add_argument("--low", action="store_true",\
		help="Remove outlying coverage")
	parser.add_argument("--plot", action='store_true',\
		help="Plot Statistics")
	args = parser.parse_args()
	if os.path.splitext(args.F)[1] in ['.fasta','.fa']:
		fai = args.F+'.fai'
	else:
		sys.exit("Please specify a fasta file\n")
	fList = parseIN(args.infile)
	#####################################################
	# Convert BAMs to Bedgraph files
	#####################################################
	makeBedgraph(fList, args.F, args.S, args.C.lower())
	#####################################################
	# Calculate replication ratio
	#####################################################
	chromDict = readFAI(fai)
	run_logFC(fList, args.norm, args.use, args.low, args.plot)
	#####################################################
	# Apply Haar Wavelet
	#####################################################
	smooth(args.L, fList, chromDict, args.use)
	gc.collect()
	#####################################################
	# Perform Segmentation
	#####################################################
	makeGFF(fList, chromDict, args.L, args.S, args.plot, args.rep, args.use, args.T, args.P, args.seg)
	print "Done"

def memAvail(p=0.8):
	mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_AVPHYS_PAGES')  # e.g. 4015976448
	mem_gib = mem_bytes/(1024.**3)  # e.g. 3.74
	return int(mem_gib*p)

def run_logFC(fList, normMethod, use, low, plotCoverage):# thresh=2.5):
	if use == 'log':
		fSuff = '_logFC.bedgraph'
	else:
		fSuff = '_ratio.bedgraph'
	beds = map(lambda y: y[1]+".bedgraph", fList)
	if not os.path.exists(fList[1][1]+fSuff):
		print "Making LogFold Files From:", beds
		L, naVals = processFiles(beds)
		gc.collect()
		#threshVals(naVals[0,:], thresh) # raise low counts to minimal coverage in G1
		# Dont think I need to do this anymore
		# Remove coverage levels below 2.5% and above 97.5%
		times = map(lambda y: y[1], fList)
		if low:
			removeLowCoverage(naVals, times, plotCoverage)
		if normMethod == 'DESeq':
			normVals = DENormalize(naVals)
		elif normMethod == 'Coverage':
			normVals = covNorm(naVals)
		else:
			sys.exit("%s is not a valid normalization method."%(normMethod))
		gtZero = normVals[0,:] > 0.0
		whereZero = normVals[0,:] == 0.0
		for bIndex in xrange(1,len(beds)):
			outName = fList[bIndex][1]+fSuff
			if not os.path.exists(outName):
				print "Making %s"%(outName)
				OF = open(outName,'w')
				if use == 'log':
					subV = (normVals[bIndex,:]+1.0)/(normVals[0,:]+1.0)
					lVals = np.log2(subV)
				else:
					lVals = np.zeros(normVals.shape[1])
					lVals[gtZero] = normVals[bIndex,gtZero]/(normVals[0,gtZero])
					# decided to remove the control=zero areas
					#lVals[whereZero] = normVals[bIndex,whereZero]
				for i in xrange(len(L[0])):
					outStr = '%s\t%i\t%i\t%.4f\n' % (L[0][i], L[1][i], L[2][i], lVals[i])
					OF.write(outStr)
				OF.close()
		if use == "log": del subV
		del lVals, L, normVals
	gc.collect()

def removeLowCoverage(vals, names, plotCoverage, cut=0.05):
	'''
	Remove  2.5% > log(coverage) > 97.5%
	
	I originally investigated the possibility of using
	a poisson distribution or a skewed normal, but the
	poisson wasn't appropriate for these coverage levels
	even with the sqrt transform, and the skewed normal
	was skewed too much by the maximum values. A regular
	normal distribution on the sqrt transform ended up
	matching the data the best.

	Normal distribution also didn't fit the data well. The
	lows were negative and the upper was too large because
	of the skew. Decded to just take the percentiles of the
	raw data. This was bad because it just took the the
	first and last 2.5% of the values.

	Found out the data looks poisson with a log transform.
	'''
	lowerP = cut/2.0
	upperP = 1.0-lowerP
	print "Removing %.1f%% > coverage > %.1f%%" % (lowerP*100, upperP*100)
	for i in range(vals.shape[0]):
		plt.figure(i)
		logNZ = np.log(vals[i,vals[i,:] > 0])
		n, b, p = plt.hist(logNZ, label="Data", bins=100, normed=True)
		yMax = np.max(n)
		gShape, gLoc, gScale = ss.gamma.fit(logNZ)
		plt.xlim(0,max(logNZ))
		plt.ylim(0,yMax)
		X = np.arange(0, max(logNZ), 0.1)
		lowerCut, upperCut = ss.gamma.ppf([lowerP, upperP], gShape, gLoc, gScale)
		if plotCoverage:
			plt.plot(X, ss.gamma.pdf(X,gShape, gLoc, gScale), label="Gamma Fit")
			plt.legend()
			plt.fill_between(X, 0, yMax, where=np.logical_or(X<lowerCut, X>upperCut), color='gray', alpha=0.6)
			plt.title("%s Coverage Cut"%(names[i]))
			plt.xlabel("log(Reads/Bin)")
			plt.ylabel("%")
			plt.savefig('%s_coverage_cut.png'%(names[i]))	
		vals[i,vals[i,:] < np.exp(lowerCut)] = 0
		vals[i,vals[i,:] > np.exp(upperCut)] = 0

def covNorm(vals):
	'''
	Normalize the counts by setting genomic coverage to 1.

	>>> covNorm(np.ones(6).reshape(2,3))
	array([[ 1.,  1.,  1.],
               [ 1.,  1.,  1.]])
	>>> covNorm(np.arange(6).reshape(2,3))
	array([[ 0.  ,  1.  ,  2.  ],
	       [ 0.75,  1.  ,  1.25]])
	'''
	sums = np.sum(vals,axis=1,dtype=np.long)
	nVals = np.float(vals.shape[1])
	return vals*(float(nVals)/sums.reshape(len(sums),1))

def DENormalize(vals):
	'''
	Normalize the counts using DESeq's method

	>>> np.round(DENormalize(np.array([[1,2,3],[4,5,6]])),2)
	array([[ 1.58,  3.16,  4.74],
	       [ 2.53,  3.16,  3.79]])
	'''
	sf = sizeFactors(vals)
	nA = np.array(vals/np.matrix(sf).T)
	return nA

def processFiles(files):
	p = Pool(cpu_count())
	vals = []
	for i in p.imap(parseVals, files):
		vals.append(i)
	p.close()
	p.join()
	del p
	gc.collect()
	locations = parseLocs(files[0])
	return (locations, np.array(vals,dtype=np.float32))

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

def smooth(level, fList, chromDict, use):
	sortedChroms = sorted(chromDict.keys()[:])
	if use == 'log':
		fSuff = '_logFC.bedgraph'
	else:
		fSuff = '_ratio.bedgraph'
	beds = map(lambda y: y[1]+fSuff, fList[1:])
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
	array([ 2.        ,  3.16227766,  4.24264069])
	>>> geometricMean(np.array([[0.1, 1, 2],[10,1,0.2]]))
	array([ 1.        ,  1.        ,  0.63245553])
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
		fList.append((tuple(tmp[1:]), tmp[0]))
	return fList

def sortWorker(inFile, chrom, nCols, conn, OF):
	bedHeap = []
	bPipe = sp.Popen('grep "^%s\s" %s'%(chrom, inFile), shell=True, stdout=sp.PIPE).stdout
	if nCols == 3:
		for line in bPipe:
			tmp = line.rstrip('\n').split('\t')
			s = int(tmp[1])
			e = int(tmp[2])
			bedHeap.append((s,e))
	else:
		for line in bPipe:
			tmp = line.rstrip('\n').split('\t')
			s = int(tmp[1])
			e = int(tmp[2])
			o = '\t'.join(tmp[3:])
			bedHeap.append((s,e,o))
	bPipe.close()
	bedHeap.sort(key=lambda y: y[0])
	msg = conn.recv()
	if nCols == 3:
		for s,e in bedHeap:
			OF.write("%s\t%i\t%i\n"%(chrom,s,e))
	else:
		for s,e,o in bedHeap:
			OF.write("%s\t%i\t%i\t%s\n"%(chrom,s,e,o))
	OF.flush()

def sortBED(inFile, sortedChroms, outFile=""):
	#sortedChroms = getSortedChroms(inFile)
	nCols = open(inFile,'r').readline().count('\t')+1
	if nCols < 3:
		sys.exit("Not a bed")
	if outFile:
		OF = open(outFile,'w')
	else:
		OF = sys.stdout
	pipes = []
	procs = []
	for chrom in sortedChroms:
		pConn, cConn = Pipe()
		pipes.append(pConn)
		procs.append(Process(target=sortWorker, args=(inFile, chrom, nCols, cConn, OF)))
		procs[-1].start()
	for i in xrange(len(sortedChroms)):
		pipes[i].send(True)
		procs[i].join()
	OF.close()

def bamtobed(bam):
	bamBed = os.path.splitext(bam)[0]+'.bed'
	if not os.path.exists(bamBed):
		os.system("bedtools bamtobed -i %s | cut -f 1-3 > %s"%(bam, bamBed))
	return True

def makeBedgraph(fList, fasta, size, replicates):
	fai = fasta+".fai"
	faiBed = "%s.%i.bed"%(fasta,size)
	sortedChroms = sorted(readFAI(fai).keys())
	memG = memAvail()
	if not os.path.exists(fai):
		print "Making fai index"
		os.system("samtools faidx %s"%(fasta))
	if not os.path.exists(faiBed):
		print "Making %ibp window bed from FAI"%(size)
		os.system("bedtools makewindows -g %s -w %i | sort -S %iG -k1,1 -k2,2n > %s"%(fai,size,memG,faiBed))
	## generate sorted beds
	bamList = []
	for bams, prefix in fList:
		for bam in bams:
			bamBed = os.path.splitext(bam)[0]+'.bed'
			if not os.path.exists(bamBed):
				bamList.append(bam)
	if len(bamList) > 0:
		p = Pool(min((cpu_count(),len(bamList))))
		ret = p.map(bamtobed, bamList)
		p.close()
		p.join()
		for bam in bamList:
			bamBed = os.path.splitext(bam)[0]+'.bed'
			print "Generating %s"%(bamBed)
			sortBED(bamBed, sortedChroms, 'tmp.bed')
			os.system("mv tmp.bed %s"%(bamBed))
	## generate bedgraphs
	for bams, prefix in fList:
		finalBG = prefix+'.bedgraph'
		if not os.path.exists(finalBG):
			if len(bams) == 1:
				bamBed = os.path.splitext(bams[0])[0]+'.bed'
				print "Generating intersect for %s"%(bamBed)
				os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%(faiBed,bamBed,finalBG))
			else:
				bgs = []
				for bam in bams:
					bamBase = os.path.splitext(bam)[0]
					bamBed = bamBase+'.bed'
					print "Generating intersect for %s"%(bamBed)
					bg = bamBase+'.bedgraph'
					os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%(faiBed,bamBed,bg))
					bgs.append(bg)
				bgStr = ' '.join(bgs)
				print "Merging\n- "+'\n- '.join(bgs)
				os.system("sort -m -S %iG -k1,1 -k2,2n %s | bedtools map -a %s -b stdin -c 4 -o %s > %s && rm %s"%(memG, bgStr, faiBed, replicates, finalBG, bgStr))

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

def fMissingCoverage(t, allSignal):
	max = allSignal.shape[1]
	ret= np.mean(np.sum(np.any(allSignal > t, axis=0)))/float(max)
	return ret

def plotCoverage(dX, d1, thresh, intF, chrom):
	plt.figure(1)
	plt.subplot(211)
	plt.plot(dX,intF(dX))
	plt.ylim(0,1)
	plt.axvline(x=thresh,color="red")
	plt.title("Cubic Interpolation of %s Coverage @ %.2f"%(chrom, thresh))
	plt.ylabel("Fraction of Chromosome")
	plt.subplot(212)
	plt.plot(dX,d1)
	plt.axvline(x=thresh,color="red")
	plt.title("Derivative of Interpolation for %s @ %.2f"%(chrom, thresh))
	plt.ylabel("Fraction of Chromosome")
	plt.xlabel("Threshold")
	plt.savefig("%s_fig.png"%(chrom))
	plt.clf()

def makeGFF(fList, chromDict, level, S, plotCov, threshMethod, use, thresh, pCut, segMeth):
	sortedChroms = sorted(chromDict.keys()[:])
	fSuff = {'log':'logFC', 'ratio':'ratio'}
	beds = map(lambda y: "%s_%s_%i.smooth.bedgraph"%(y[1], fSuff[use], level), fList[1:])
	names = map(lambda y: y[1], fList[1:])
	for name in names:
		outGFF = '%s_%s_%i.smooth.gff3'%(name, fSuff[use], level)
		open(outGFF,'w').write('##gff-version 3\n#track name="%s LogFold %ibp" gffTags=on\n' % (name,S))
		open(fSuff[use]+"_segmentation.gff3",'w').write('##gff-version 3\n#track name="Segmentation %ibp" gffTags=on\n'%(S))
	print "Parsing :",beds
	L, vals = processFiles(beds)
	locDict = parseLocations(L[0])
	for chrom in sortedChroms:
		s,e = locDict[chrom]
		tmpChr = L[0][s:e]
		tmpS = L[1][s:e]
		tmpE = L[2][s:e]
		maskM = np.zeros((len(names),len(tmpChr)),dtype=np.bool)
		allSignal = vals[:len(beds),s:e]
		if threshMethod == "threshold":
			pass
		elif threshMethod == "auto":
			asMin = np.min(allSignal)
			asMax = np.max(allSignal)
			X = np.arange(asMin-0.5, asMax+0.5, 0.1)
			vFunc = np.vectorize(fMissingCoverage, excluded=['allSignal'])
			Y = vFunc(X,allSignal=allSignal)
			intF = scipy.interpolate.interp1d(X,Y,kind="cubic")
			dX = np.arange(asMin,asMax,0.01)
			d1 = scipy.misc.derivative(intF,dX, dx=0.05,n=1)
			try:
				locMax = np.argmax(np.abs(d1))
				lowerLocs = np.where(np.abs(d1) < 0.1)[0]
				belowLocal = lowerLocs[lowerLocs < locMax]
				thresh = dX[np.max(belowLocal)]
			except:
				thresh = 0.0
			print "%s replication threshold: %.2f"%(chrom, thresh)
			if plotCov:
				plotCoverage(dX, d1, thresh, intF, chrom)
		elif threshMethod == "percent":
			thresh = np.percentile(allSignal, pCut)
			#bThresh = np.percentile(allSignal, pCut, axis=1) 
		else:
			sys.exit("invalid threshold method")
		for i in xrange(len(beds)):
			outGFF = '%s_%s_%i.smooth.gff3'%(names[i], fSuff[use], level)
			outSignal = vals[i,s:e]
			bMask = np.zeros(len(outSignal), dtype=np.bool)
			bMask[outSignal > thresh] = 1
			maskM[i,:]=bMask[:]
			rBounds = calcRegionBounds(bMask)
			OF = open(outGFF,'a')
			count = 1
			for rS,rE in rBounds:
				OF.write("%s\t.\tgene\t%i\t%i\t.\t.\t.\tID=gene%i;color=%s;\n" % (tmpChr[rS], tmpS[rS]+1, tmpE[rE], count, myColors[i]))
				count += 1
			OF.close()
		## Segmentation
		if segMeth == "proportion":
			for i in xrange(np.shape(maskM)[1]):
				maskRow = maskM[:,i]
				vRow = np.copy(allSignal[:,i])
				vRow[maskRow == 0] = 0.0
				rowSum = np.sum(maskRow)
				if rowSum > 1:
					maskM[:,i] = classProportion(vRow)
		elif segMeth == "binary":
			pass
		else:
			sys.exit("invalid segmentation method: "+segMeth)
		OS = open(fSuff[use]+'_segmentation.gff3','a')
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

def classProportion(dataRow, emlSize = 0.5):
	'''
	Calculates the segmentation class based on the
	proportion of replication.

	>>> classProportion(np.array([1.0,0.49,0.51]))
	array([ True, False,  True], dtype=bool)
	'''
	maxVal = np.float(np.max(dataRow))
	maxInd = np.argmax(dataRow)
	tRow = dataRow/float(maxVal)
	nd = len(dataRow)
	ndm1 = nd-1
	emlVol = emlSize**ndm1
	emlCut = 1.0-emlSize
	otherSpaces = 2.0**ndm1-1.0
	firstCut = ((1.0-emlVol)/otherSpaces)**(1.0/ndm1)
	out = np.zeros(nd, dtype=np.bool)
	if np.all(tRow > emlCut):
		return np.ones(nd, dtype=np.bool)
	out[tRow > emlCut] = 1
	for i in range(ndm1):
		for j in range(i+1,nd):
			if np.all(out[[i,j]]):
				pass
			elif np.all(tRow[[i,j]] > firstCut):
				if tRow[i] > tRow[j]:
					out[i] = 1
				elif tRow[j] > tRow[i]:
					out[j] = 1
			else:
				if tRow[i] > firstCut:
					out[i] = 1
				if tRow[j] > firstCut:
					out[j] = 1
	return out

def powerSet(L):
	'''
	Removes the empty list from powerset

	>>> powerSet([0,1])
	[[1], [0], [1, 0]]
	'''
	def powerHelp(A):
		'''
		Builds powerset recursively.

		>>> powerHelp([0,1])
		[[], [1], [0], [1, 0]]
		'''
		if not A:
			return [[]]
		ret = powerHelp(A[1:])
		return ret+[i+[A[0]] for i in ret]
	return powerHelp(L)[1:]

def haarVec(vals, p=90):
	'''
	Performs a Haar wavelet transformation on the array,
	removes the lower p% of the values, and then
	inverse transforms the remaining values.

	>>> haarVec([1,2,5,1,2,3,4,0,7,3,2,1])
	array([ 2.25,  2.25,  2.25,  2.25,  2.25,  2.25,  2.25,  2.25,  5.  ,
	        5.  ,  1.5 ,  1.5 ])
	'''
	import pywt
	hC = pywt.wavedec(vals,'haar')
	## Updated to not count the average (val 0)
	cutVal = np.percentile(np.abs(np.concatenate(hC[1:])), p)
	for A in hC[1:]: ## skip 0 (overall average)
		A[np.abs(A) <= cutVal] = 0
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
