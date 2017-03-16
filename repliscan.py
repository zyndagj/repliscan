#!/usr/bin/env python

import os, sys, math, argparse, gc, colorsys
import numpy as np
from multiprocessing import Pool, cpu_count, Process, Pipe
import subprocess as sp
from glob import glob
from copy import deepcopy
from array import array
from itertools import compress
import scipy.optimize
import scipy.interpolate
import scipy.misc
import scipy.stats as ss
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
	parser.add_argument("config", metavar="FILE", help="File with list of bams")
	parser.add_argument("-r", "--ref", metavar='FASTA',help="Reference fasta", required=True)
	parser.add_argument("-l", "--level", metavar='INT', help=\
		"Haar smoothing level (Default: %(default)s)", default=3, type=int)
	parser.add_argument("-w", "--window", metavar='INT', help=\
		"Analysis bin size in base pairs (Default: %(default)s)", default=1000, type=int)
	parser.add_argument("-a", "--aggregate", metavar='STR', help=\
		"Replicate agregation method (sum|median|mean|min|max) (Default: %(default)s)", default="sum",\
		type=argChecker(('sum','median','mean','min','max'), 'reducer').check)
	#parser.add_argument("-n", "--norm", metavar='STR', help=\
	#	"Normalization Method (DESeq|Coverage) (Default: %(default)s)", default="Coverage",\
	#	type=argChecker(('DESeq','Coverage'),'normalization method').check)
	parser.add_argument("-t", "--threshold", metavar='STR', help=\
		"Replication threshold method (value|auto|percent) (Default: %(default)s)", default="auto",\
		type=argChecker(('value','auto','percent'),'replication method').check)
	parser.add_argument("-S", "--scope", metavar='STR', help=\
		"Replication scope (chromosome|genome) (Default: %(default)s)", default="chromosome",\
		type=argChecker(('chromosome','genome'),'replication scope').check)
	parser.add_argument("-v","--value", metavar='FLOAT', help=\
		"Explicit replication threshold value of [%(default)s] when using '-t value'", default=1.0, type=float)
	parser.add_argument("--prep", metavar='FLOAT', help=\
		"Consider the lowest [%(default)s]%% of the signal to be noise and not replicating when using '-t percent'", default=2.0, type=float)
	parser.add_argument("-c", "--classifier", metavar='STR', help=\
		"Segmentation classification method (binary|proportion) (Default: %(default)s)", default=\
		"proportion", type=argChecker(('binary','proportion'),'segmentation method').check)
	parser.add_argument("-R", "--remove", metavar="STR", help=\
		"Outlying data to remove (none|lognGamma|sqrtGamma|norm|whiskers|percentile) (Default: %(default)s)", default=\
		"norm", type=argChecker(("none","lognGamma","sqrtGamma","norm","whiskers","percentile"), "outlying method").check)
	parser.add_argument("--pcut", metavar='FLOAT', help=\
		"Remove the upper and lower [%(default)s]%% of the data when using '-R percentile'", default=2.5, type=float)
	parser.add_argument("--log", action="store_true", help=\
		"Apply log transform to sequenceability ratio (Default: False)")
	parser.add_argument("-f","--force", action="store_true",\
		help="Force the re-generation of all files")
	parser.add_argument("--plot", action='store_true',\
		help="Plot Statistics")
	args = parser.parse_args()
	if os.path.splitext(args.ref)[1] in ['.fasta','.fa']:
		fai = args.ref+'.fai'
	else:
		sys.exit("Please specify a fasta file\n")
	# Set force variable
	global force
	force = args.force
	#####################################################
	# Parse exp config
	#####################################################
	fList = parseIN(args.config)
	#####################################################
	# Convert BAMs to Bedgraph files
	#####################################################
	#L, normVals = makeBedgraph(fList, args.ref, args.window, args.aggregate.lower(), args.norm, args.remove, args.plot)
	L, normVals = makeBedgraph(fList, args.ref, args.window, args.aggregate.lower(), 'Coverage', args.remove, args.plot, args.pcut)
	#####################################################
	# Calculate replication ratio
	#####################################################
	run_logFC(fList, args.log, L, normVals)
	#####################################################
	# Apply Haar Wavelet
	#####################################################
	chromDict = readFAI(fai, args.window)
	smooth(args.level, fList, chromDict, args.log)
	gc.collect()
	#####################################################
	# Perform Segmentation
	#####################################################
	makeGFF(fList, chromDict, args.level, args.window, args.plot, args.threshold, args.scope, args.log, args.value, args.prep, args.classifier)
	print "Done"

def memAvail(p=0.8):
	mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_AVPHYS_PAGES')  # e.g. 4015976448
	mem_gib = mem_bytes/(1024.**3)  # e.g. 3.74
	return int(mem_gib*p)

def normalizer(method, vals):
	'''
	Normalization superfunction.

	>>> normalizer('Coverage',np.ones(6).reshape(2,3))
	Normalizing signales using Coverage
	array([[ 1.,  1.,  1.],
	       [ 1.,  1.,  1.]])
	>>> np.round(normalizer('DESeq',np.array([[1,2,3],[4,5,6]])),2)
	Normalizing signales using DESeq
	array([[ 1.58,  3.16,  4.74],
	       [ 2.53,  3.16,  3.79]])
	'''
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
		sums = np.sum(vals,axis=1)
		#sums = np.sum(vals,axis=1, dtype=np.long)
		nVals = np.float(vals.shape[1])
		return vals*(nVals/sums.reshape(len(sums),1))
	def DENormalize(vals):
		'''
		Normalize the counts using DESeq's method

		>>> np.round(DENormalize(np.array([[1,2,3],[4,5,6]])),2)
		array([[ 1.58,  3.16,  4.74],
		       [ 2.53,  3.16,  3.79]])
		'''
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
		sf = sizeFactors(vals)
		nA = np.array(vals/np.matrix(sf).T)
		return nA
	D = {'DESeq':DENormalize, 'Coverage':covNorm}
	if method in D:
		print "Normalizing signales using %s"%(method)
		return D[method](vals)
	else:
		sys.exit("%s is not a valid normalization method."%(normMethod))

def run_logFC(fList, useLog, L, normVals):# thresh=2.5):
	if useLog:
		fSuff = '_logFC.bedgraph'
	else:
		fSuff = '_ratio.bedgraph'
	beds = map(lambda y: y[1]+"_norm.bedgraph", fList)
	print "Making LogFold Files From:", beds
	#L, normVals = processFiles(beds)
	#threshVals(naVals[0,:], thresh) # raise low counts to minimal coverage in G1
	# Dont think I need to do this anymore
	# Remove coverage levels below 2.5% and above 97.5%
	gtZero = normVals[0,:] > 0.0
	for bIndex in xrange(1,len(beds)):
		outName = fList[bIndex][1]+fSuff
		if not os.path.exists(outName) or force:
			print "Making %s"%(outName)
			OF = open(outName,'w')
			if useLog:
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

def removeOutlying(vals, names, plotCoverage, method, pOut):
	'''
	Removes outlying coverage values from the vals matrix in-place.
	'''
	if method == 'none': return
	func = {'norm':calcNorm, 'sqrtGamma':sqrtGamma, 'whiskers':calcWhisk, 'lognGamma':lognGamma, 'percentile':calcPercentile}
	N = len(names)
	p = Pool(min((cpu_count(), N)))
	cuts = p.map(func[method], zip(vals, names, [plotCoverage]*N, [pOut]*N))
	p.close()
	p.join()
	for i in range(N):
		vals[i,vals[i,:] < cuts[i][0]] = 0
		vals[i,vals[i,:] > cuts[i][1]] = 0
def lognGamma(argList):
	V, name, plotCoverage, pOut = argList
	cut = 0.05
	lowerP, upperP = cut/2.0, 1.0-cut/2.0
	transV = np.log(V[V > 0])
	fits = ss.gamma.fit(transV)
	lowerCut, upperCut = ss.gamma.ppf([lowerP, upperP], *fits)
	if plotCoverage:
		plt.figure(figsize=(10,3))
		X = plotCut(transV, name, lowerCut, upperCut)
		plt.plot(X, ss.gamma.pdf(X, *fits), label="Gamma Fit", color="#FFC107", lw=2)
		plt.legend()
		plt.xlabel("log(Reads/Bin)")
		plt.tight_layout()
		plt.savefig('%s_coverage_cut.png'%(name))
	return np.exp((lowerCut, upperCut))
def plotCut(V, name, lowerCut, upperCut):
	n, b, p = plt.hist(V, label="Data", bins=100, normed=True, color="#03A9F4")
	minV, maxV, maxY = min(V), max(V), np.max(n)
	plt.xlim(minV,maxV)
	plt.ylim(0,maxY)
	X = np.arange(minV, maxV+0.1, 0.1)
	plt.fill_between(X, 0, maxY, where=np.logical_or(X<lowerCut, X>upperCut), color='gray', alpha=0.6)
	plt.title("%s Coverage Cut"%(name))
	plt.ylabel("Fraction of Reads")
	return X
def calcWhisk(argList):
	V, name, plotCoverage, pOut = argList
	transV = np.log(V[V > 0])
	q1, q3 = np.percentile(transV, (25,75))
	iqr = q3-q1
	lWhisk = q1-1.5*iqr
	uWhisk = q3+1.5*iqr
	if plotCoverage:
		plt.figure(figsize=(10,3))
		X = plotCut(transV, name, lWhisk, uWhisk)
		plt.xlabel("log(Reads/Bin)")
		plt.tight_layout()
		plt.savefig('%s_coverage_cut.png'%(name))
	return np.exp((lWhisk, uWhisk))
def calcPercentile(argList):
	V, name, plotCoverage, pOut = argList
	transV = np.log(V[V > 0])
	q1, q3 = np.percentile(transV, (pOut, 1-pOut))
	if plotCoverage:
		plt.figure(figsize=(10,3))
		X = plotCut(transV, name, q1, q3)
		plt.xlabel("log(Reads/Bin)")
		plt.tight_layout()
		plt.savefig('%s_coverage_cut.png'%(name))
	return np.exp((q1, q3))
def sqrtGamma(argList):
	V, name, plotCoverage, pOut = argList
	transV = np.sqrt(V)
	q1, q3 = np.percentile(transV, (25,75))
	iqr = q3-q1
	lWhisk = q1-1.5*iqr
	uWhisk = q3+1.5*iqr
	lP, uP = 0.025, 0.975
	transV = transV[transV < uWhisk]
	transV = transV[transV > lWhisk]
	fits = ss.gamma.fit(transV)
	lowerCut, upperCut = ss.gamma.ppf([lP, uP], *fits)
	if plotCoverage:
		plt.figure(figsize=(10,3))
		X = plotCut(transV, name, lowerCut, upperCut)
		plt.plot(X, ss.gamma.pdf(X, *fits), label="Gamma Fit", color="#FFC107", lw=2)
		plt.legend()
		plt.xlabel("sqrt(Reads/Bin)")
		plt.tight_layout()
		plt.savefig('%s_coverage_cut.png'%(name))
	return np.power((lowerCut, upperCut), 2)
def calcNorm(argList):
	V, name, plotCoverage, pOut = argList
	transV = np.log(V[V > 0])
	q1, q3 = np.percentile(transV, (25,75))
	lP, uP = 0.025, 0.975
	fits = ss.norm.fit(transV)
	lowerCut, upperCut = ss.norm.ppf([lP, uP], *fits)
	if plotCoverage:
		plt.figure(figsize=(10,3))
		X = plotCut(transV, name, lowerCut, upperCut)
		plt.plot(X, ss.norm.pdf(X, *fits), label="Normal Fit", color="#FFC107", lw=2)
		plt.legend()
		plt.xlabel("log(Reads/Bin)")
		plt.tight_layout()
		plt.savefig('%s_coverage_cut.png'%(name))
	return np.exp((lowerCut, upperCut))

def processFiles(files):
	'''
	Returs the locations and values for a collection of bedgraph files.
	'''
	vals = runPool(parseVals, files)
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

def smooth(level, fList, chromDict, useLog):
	sortedChroms = sorted(chromDict.keys()[:])
	if useLog:
		fSuff = '_logFC.bedgraph'
	else:
		fSuff = '_ratio.bedgraph'
	beds = map(lambda y: y[1]+fSuff, fList[1:])
	print "Smoothing:",beds
	for bed in beds:
		outFile = '%s_%i.smooth.bedgraph'%(bed.rstrip('.bedgraph'),level)
		if not os.path.exists(outFile) or force:
			os.system('[ -e %s ] && rm %s'%(outFile, outFile))
			for chr in sortedChroms:
				locCmd = "grep '^%s\s' %s | cut -f 1-3"%(chr, bed)
				valCmd = "grep '^%s\s' %s | cut -f 4 | wavelets --level %i --to-stdout --boundary reflected --filter Haar -" % (chr, bed, level)
<<<<<<< HEAD
				os.system('bash -c "paste <( %s ) <( %s ) >> %s"'%(locCmd, valCmd, outFile))
=======
				awkCmd = "awk '{if (NF == 4) print \$0;}'"
				os.system('bash -c "paste <( %s ) <( %s ) | %s >> %s"'%(locCmd, valCmd, awkCmd, outFile))
>>>>>>> b7feb8a0ec3c661af78dd1af9bdf8c6433086082
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

def readFAI(inFile, windowSize):
	'''
	Returns (length, seek, read length, readlength + newline)
	'''
	#Pt     140384  50      60      61
	chromDict = {}
	for line in open(inFile,'r'):
		tmp = line.split('\t')
		#chromDict[tmp[0]] = tuple(map(int, tmp[1:]))
		if int(tmp[1])/windowSize > 30:
			chromDict[tmp[0]] = int(tmp[1])
	return chromDict

def parseIN(inFile):
	fList = []
	for line in open(inFile,'r'):
		tmp = line.rstrip('\n').split('\t')
		fList.append((tuple(filter(lambda x: x != '', tmp[1:])), tmp[0]))
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
	if not os.path.exists(bamBed) or force:
		os.system("export LC_ALL=C; bedtools bamtobed -i %s | cut -f 1-3 | sort -S 1G -k1,1 -k2,2n > %s"%(bam, bamBed))
	return bamBed

def bedtobedgraph(argList):
	faiBed, bed = argList
	bedgraph = os.path.splitext(bed)[0]+'.bedgraph'
	if not os.path.exists(bedgraph) or force:
		os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%(faiBed, bed, bedgraph))
	return bedgraph

def runPool(func, args):
	'''
	Maps a function to arguments by spawning and joining processor pool.
	'''
	p = Pool(min((cpu_count(),len(args))))
	ret = p.map(func, args)
	p.close()
	p.join()
	return ret

def makeBedgraph(fList, fasta, size, aggMethod, normMethod, removeWhat, plotCoverage, pOut):
	def normBGs(files, normMethod):
		L, vals = processFiles(files)
		normVals = normalizer(normMethod, vals)
		normList = map(lambda x: os.path.splitext(x)[0]+'_norm.bedgraph', files)
		writeVals(L, normVals, normList)
		return normList
	fai = fasta+".fai"
	faiBed = "%s.%i.bed"%(fasta,size)
	memG = memAvail()
	##########################################
	# Make reference windows
	##########################################
	if not os.path.exists(fai) or force:
		print "Making fai index"
		os.system("samtools faidx %s"%(fasta))
	if not os.path.exists(faiBed) or force:
		print "Making %ibp window bed from FAI"%(size)
		windowCMD="export LC_ALL=C; bedtools makewindows -g %s -w %i | sort -S %iG -k1,1 -k2,2n > %s"%(fai,size,memG,faiBed)
		os.system(windowCMD)
	##########################################
	# Convert bam to sorted bed
	##########################################
	print "Converting bam to bed"
	bamList = [bam for bams, prefix in fList for bam in bams]
	bedList = runPool(bamtobed, bamList)
	##########################################
	# Convert bed to bedgraph
	##########################################
	print "Converting bed to bedgraph"
	bgList = runPool(bedtobedgraph, [(faiBed, bed) for bed in bedList])
	##########################################
	# Normalize replicates unless using sum
	##########################################
	if aggMethod == "sum":
		normBG = bgList
	else:
		normBG = normBGs(bgList, normMethod)
	##########################################
	# Aggregate replicates
	##########################################
	print "Aggregating replicates"
	prefixList = runPool(agregate, [(bams, prefix, faiBed, aggMethod) for bams, prefix in fList])
	L, vals = processFiles(prefixList)
	##########################################
	# Remove outlying coverage
	##########################################
	times = map(lambda y: y[1], fList)
	removeOutlying(vals, times, plotCoverage, removeWhat, pOut)
	##########################################
	# Normalized agregated files
	##########################################
	normVals = normalizer(normMethod, vals)
	normPrefix = map(lambda x: os.path.splitext(x)[0]+'_norm.bedgraph', prefixList)
	writeVals(L, normVals, normPrefix)
	return (L, normVals)

def writeVals(L, vals, files):
	'''
	Iterates of a list values and writes each array to a bedgraph file

	Parameters
	------------------------
	L	(3,x) tuple of bedgraph coordinates
	vals	(N,x) numpy array of values
	files	(N) tuple of output file names

	Returns
	------------------------
	'''
	c,s,l = L
	for i in range(len(files)):
		if not os.path.exists(files[i]) or force:
			if len(vals[i]) != len(l): sys.exit("Lengths don't match")
			open(files[i], 'w').write('\n'.join(map(lambda w,x,y,z: '\t'.join(map(str, (w,x,y,z))), c,s,l,vals[i]))+'\n')

def agregate(argList):
	'''
	Agregates bedgraph replicates using bedtools map.

	Parameters
	------------------------
	argList		Tuple of (bams, prefix, faiBed, method)

	Returns
	------------------------
	finalBG		Name of the output bedgraph
	'''
	bams, prefix, faiBed, method = argList
	finalBG = '%s.bedgraph'%(prefix)
	if method == 'sum':
		bgs = ['%s.bedgraph'%(os.path.splitext(bam)[0]) for bam in bams]
	else:
		bgs = ['%s_norm.bedgraph'%(os.path.splitext(bam)[0]) for bam in bams]
	if not os.path.exists(finalBG) or force:
		if len(bgs) == 1:
			os.system('ln -fs %s %s'%(bgs[0], finalBG))
		else:
			print "Merging with %s\n- "%(method)+'\n- '.join(bgs)
			bgStr = ' '.join(bgs)
			mapStr = "export LC_ALL=C; sort -m -S 1G -k1,1 -k2,2n %s | bedtools map -a %s -b - -c 4 -o %s > %s"%(bgStr, faiBed, method, finalBG)
			os.system(mapStr)
	return finalBG

def parseLocations(chroms):
	'''
	Returns a dictionary of the start and end indexes of each chromosome that correspond to the location arrays.

	Parameters
	------------------------
	chroms		List of chromosomes from L

	Returns
	------------------------
	locDict		Dictionary of {chrom:(start,end)} inidces
	'''
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

def calcThreshold(vals, locDict, threshMethod, scope, T, pCut, plotCov):
	'''
	>>> calcThreshold([0,2,1], {1:(2,3), 2:(4,5)}, 'value', 'genome', 1.0, 2.5, True)
	{1: 1.0, 2: 1.0}
	'''
	if threshMethod == 'value':
		return dict.fromkeys(locDict, T)
	if scope == 'genome':
		if threshMethod == 'auto':
			thresh = autoThresh(vals, 'genome', plotCov)
		elif threshMethod == 'percent':
			thresh = perThresh(vals, pCut)
		else:
			sys.exit("\nBad thresh method for genome coverage: %s\n"%(threshMethod))
		return dict.fromkeys(locDict, thresh)
	elif scope == 'chromosome':
		threshDict = {}
		for chrom, (s,e) in locDict.iteritems():
			chrVals = vals[:,s:e]
			if threshMethod == 'auto':
				threshDict[chrom] = autoThresh(chrVals, chrom, plotCov)
			elif threshMethod == 'percent':
				threshDict[chrom] = perThresh(chrVals, pCut)
			else:
				sys.exit("\nBad thresh method for chromosome coverage: %s\n"%(threshMethod))
		nonNegMed = np.median(filter(lambda x: x != -1, threshDict.values()))
		for k,v in threshDict.iteritems():
			if v == -1:
				threshDict[k] = nonNegMed
		return threshDict
	else:
		sys.exit("\nBad scope: %s\n"%(scope))

def autoThresh(vals, name, plotCov):
	def fMissingCoverage(thresh, allSignal):
		'''
		float thresh
		float allSignal[time, location]

		>>> fMissingCoverage(9, np.arange(12).reshape((3,4)))
		0.5
		'''
		numBins = allSignal.shape[1]
		allGTthresh = np.any(allSignal > thresh, axis=0)
		totalGTthresh = np.sum(allGTthresh)
		return totalGTthresh/np.float(numBins)
	minVals = np.min(vals)
	maxVals = np.max(vals)
	if minVals == 0 and maxVals == 0: return 1.0
	X = np.arange(minVals-0.5, maxVals+0.5, 0.1)
	vFunc = np.vectorize(fMissingCoverage, excluded=['allSignal'])
	Y = vFunc(X,allSignal=vals)
	# Perform cubic interpolation to remove noisy local maxima
	intF = scipy.interpolate.interp1d(X,Y,kind="cubic")
	dX = np.arange(minVals,maxVals,0.01)
	# Calculate the derivative of the coverage at different thresholds
	d1 = scipy.misc.derivative(intF,dX, dx=0.05,n=1)
	dMin = np.argmin(d1)
	cutShoulder = np.where(d1 > -0.1)[0]
	try:
		belowMin = cutShoulder[cutShoulder < dMin]
		thresh = dX[np.max(belowMin)]
	except:
		#thresh = np.percentile(vals, 25)
		thresh = -1
	if plotCov:
		plotCoverage(dX, d1, thresh, intF, name)
	return thresh

def perThresh(vals, pCut):
	'''
	>>> perThresh(range(10), 40) == 3.6
	True
	'''
	return np.percentile(vals, pCut)

def int2name(val, nameList):
	'''
	>>> NL = ['S1', 'S2']
	>>> [ int2name(i, NL) for i in range(1,4) ]
	['S2', 'S1', 'S1S2']
	>>> int2name(4, NL)
	Traceback (most recent call last):
	 ...
	ValueError: Val exceeds the namelist
	>>> int2name(0, NL)
	Traceback (most recent call last):
	 ...
	ValueError: No name
	'''
	if not val:
		raise ValueError("No name")
	if val > (2**len(nameList))-1:
		raise ValueError("Val exceeds the namelist")
	binRep = map(int, np.binary_repr(val,len(nameList)))
	return ''.join(compress(nameList, binRep))

def plotVars(nameList):
	'''
	Defines the custom colors to be used for generating segmentation gffs.

	>>> plotVars(['ES','MS','LS'])
	['#FB0018', '#1A8A12', '#FFFD33', '#2250F1', '#EA3CF2', '#28C5CC', '#FAB427']
	>>> plotVars(['E','M','L'])
	['#FB0018', '#1A8A12', '#FFFD33', '#2250F1', '#EA3CF2', '#28C5CC', '#FAB427']
	>>> plotVars(['A','B','C'])
	[u'#8000ff', u'#2c7ef7', u'#2adddd', u'#80ffb4', u'#d4dd80', u'#ff7e41', u'#ff0000']
	'''
	if nameList == ['ES','MS','LS'] or nameList == ['E','M','L']:
		#print "Using ES, MS, LS colors"
		myColors = ["#FB0018","#1A8A12","#FFFD33","#2250F1","#EA3CF2","#28C5CC","#FAB427"]
		#          [ES,       MS,       LS,       MSLS,     ESLS,     ESMS,     ESMSLS]
	else:
		#print "Using rainbow colors"
		myColors = map(lambda x: matplotlib.colors.rgb2hex(x[:3]), cm.rainbow(np.linspace(0,1,2**len(nameList)-1)))
	return myColors

def makeGFF(fList, chromDict, level, S, plotCov, threshMethod, scope, useLog, thresh, pCut, segMeth):
	fSuff = {True:'logFC', False:'ratio'}[useLog]
	beds = map(lambda y: "%s_%s_%i.smooth.bedgraph"%(y[1], fSuff, level), fList[1:])
	names = map(lambda y: y[1], fList[1:])
	plotColors = plotVars(names)
	for name in names:
		outGFF = '%s_%s_%i.smooth.gff3'%(name, fSuff, level)
		open(outGFF,'w').write('##gff-version 3\n#track name="%s %ibp" gffTags=on\n' % (name,S))
	#open(fSuff[use]+"_segmentation.gff3",'w').write('##gff-version 3\n#track name="Segmentation %ibp" gffTags=on\n'%(S))
	print "Parsing :",beds
	L, vals = processFiles(beds)
	locDict = parseLocations(L[0])
	sortedChroms = sorted(locDict.keys()[:])
	print "Calculating thresholds..."
	thresholdDict = calcThreshold(vals, locDict, threshMethod, scope, thresh, pCut, plotCov)
	OS = open(fSuff+'_segmentation.gff3','w')
	OS.write('##gff-version 3\n#track name="Segmentation %ibp" gffTags=on\n'%(S))
	segCount = 1
	counts = np.ones(len(beds))
	for chrom in sortedChroms:
		s,e = locDict[chrom]
		tmpChr = L[0][s:e]
		tmpS = L[1][s:e]
		tmpE = L[2][s:e]
		maskM = np.zeros((len(names),len(tmpChr)),dtype=np.bool)
		allSignal = vals[:len(beds),s:e]
		thresh = thresholdDict[chrom]
		print "Using a threshold of %.2f for chromosome %s"%(thresh, chrom)
		for i in xrange(len(beds)):
			colorIndex = int(''.join([ '1' if i == j else '0' for j in range(len(beds))]),2)
			outGFF = '%s_%s_%i.smooth.gff3'%(names[i], fSuff, level)
			outSignal = vals[i,s:e]
			bMask = np.zeros(len(outSignal), dtype=np.bool)
			bMask[outSignal > thresh] = 1
			maskM[i,:]=bMask[:]
			rBounds = calcRegionBounds(bMask)
			OF = open(outGFF,'a')
			for rS,rE in rBounds:
				OF.write("%s\t.\tgene\t%i\t%i\t.\t.\t.\tID=gene%i;color=%s;\n" % (tmpChr[rS], tmpS[rS]+1, tmpE[rE], counts[i], plotColors[colorIndex-1]))
				counts += 1
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
		# Convert each binary category into a base-10 integer
		calls = np.array(map(lambda x: sum(2**np.where(x[::-1])[0]), maskM.T))
		# Calculate regions of repeting categories (including zeros)
		regions, labels = calcRegionCategories(calls)
		for (s, e), l in zip(regions, labels):
			if tmpChr[s] != chrom:
				sys.exit("Chromosomes didn't match")
			if l:
				name = int2name(l, names)
				OS.write("%s\t.\tgene\t%i\t%s\t.\t.\t.\tID=gene%i;Name=%s;color=%s;\n" % (chrom, tmpS[s]+1, tmpE[e], segCount, name, plotColors[l-1]))
				segCount += 1
	OS.close()

def classProportion(dataRow, emlSize = 0.5):
	'''
	Calculates the segmentation class based on the
	proportion of replication.
	
	Input
	=======================
	dataRow		numpy array of proportional values
	emlSize		size of EML region for unequal cuts
	
	Output
	=======================
	out		binary numpy array of most significant values

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

def calcRegionCategories(counter):
	'''
	Returns the new lower and upper bounds over overlapped regions, along
	with the category.

	Parameters
	=============================
	counter		Binary counter array
	
	Output
	=============================
	regions		NP array
	labels		NP array

	>>> calcRegionCategories(np.array([1,1,0,0,1,1,1,0,0]))
	(array([[0, 1],
	       [2, 3],
	       [4, 6],
	       [7, 8]]), array([1, 0, 1, 0]))
	'''
	regions = []
	labels = []
	sI = 0
	val = counter[0]
	eI = 0
	for index in range(1,len(counter)):
		if counter[index] == val:
			eI = index
		else:
			regions.append((sI,eI))
			labels.append(val)
			val = counter[index]
			sI = index
			eI = index
	regions.append((sI,eI))
	labels.append(val)
	return np.array(regions), np.array(labels)

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
