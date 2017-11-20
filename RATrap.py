#!/usr/bin/env python

import numpy as np
import argparse, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from itertools import izip_longest, compress
from operator import itemgetter

def main():
	parser = argparse.ArgumentParser(description="Finds the timing differences between two segmentation profiles.")
	parser.add_argument("-d",metavar="FLOAT", help="Minimum distance to be RAT (Default: %(default)s)", default=0.5, type=float)
	parser.add_argument("-S",metavar="INT", help="Tile Size (Default: %(default)s)", default=1000, type=int)
	parser.add_argument("-A",metavar="GFF3", help="First Segmentation Profile (mitotic)", required=True)
	parser.add_argument("-B",metavar="GFF3", help="Second Segmentation Profile (endo)", required=True)
	parser.add_argument("-T",metavar="STR", help="Times (Default: %(default)s)", default="ES,MS,LS", type=str)
	parser.add_argument("-F",metavar="FASTA", help="Reference", required=True)
	parser.add_argument("-O",metavar="BEDG", help="Output to bedgraph file", default=sys.stdout)
	parser.add_argument("--stats", action="store_true", help="Generate stats and figures")
	parser.add_argument("--fig",metavar="EXT", help="Figure type (Default: %(default)s)", default="pdf", type=str)
	parser.add_argument("--diff", action="store_true", help="Print fraction different and exit")
	args = parser.parse_args()
	if os.path.splitext(args.F)[1] in ['.fasta','.fa']:
		fai = args.F+'.fai'
		if not os.path.exists(fai): sys.exit("Please generate an faidx for %s\n"%(args.F))
	else:
		sys.exit("Please specify a fasta file\n")
	if args.stats and not args.O:
		sys.exit("Please specify an output bedgraph so stats are not hidden")
	chromDict = readFAI(fai)
	plotVars(args.T)
	sortedChroms = sorted(chromDict.keys())
	genomeA = processGenome(chromDict, args.S, args.A, args.stats, args.fig)
	genomeB = processGenome(chromDict, args.S, args.B, args.stats, args.fig)
	T, TD = (0.0, 0.0)
	for chrom in sortedChroms:
		pd, nd, n = proportionDifferent(genomeA[chrom], genomeB[chrom], args.d)
		T += n
		TD += nd
		print "%s\t%.5f"%(chrom, pd)
	print "Genome:\t%.5f"%(TD/T)
	if args.diff: return 0
	OF = open(args.O,'w',1000000)
	OF.write('#Chromosome\tstart\tend\tdistance\tA\tB\tindex-m_A\tindex-m_B\n')
	X = [[] for i in xrange(len(sortedChroms))]
	if args.stats:
		for record in compareGenomes(genomeA, genomeB, chromDict, args.d, args.S, args.stats, args.fig):
			tmp = record.split('\t')
			size = int(tmp[2])-int(tmp[1])
			X[sortedChroms.index(tmp[0])].append(size)
			OF.write(record+'\n')
		plt.figure()
		plt.boxplot(X, labels=sortedChroms, showfliers=False)
		plt.ylabel("RAT size (bp)")
		plt.xlabel("Chromosome")
		plt.title("Size Distributions of RATs")
		plt.savefig("RAT_size."+args.fig)
		plt.close()
	else:
		for record in compareGenomes(genomeA, genomeB, chromDict, args.d, args.S, args.stats, args.fig):
			OF.write(record+'\n')
	OF.close()

def helper(args):
	a,b,diff_thresh = args
	return abs(dist(a,b)) >= diff_thresh
def proportionDifferent(A, B, diff_thresh):
	'''
	>>> proportionDifferent([1,1,1],[1,1,1])
	(0.0, 0, 3)
	>>> proportionDifferent([1,1,0,0],[1,1,1,1])
	(0.5, 2, 4)
	>>> proportionDifferent([1,1,0,0],[1,1,1])
	Traceback (most recent call last):
		...
	AssertionError
	'''
	assert(len(A) == len(B))
	from multiprocessing import Pool
	from itertools import izip, repeat
	p = Pool(processes=4)
	numDifferent = sum(p.imap_unordered(helper, izip(A,B,repeat(diff_thresh,len(A))), 1000))
	p.close()
	p.join()
	return float(numDifferent)/len(A), numDifferent, len(A)

def compareGenomes(A, B, chromDict, minD, tileSize, statsFlag, figExt):
	sortedChroms = sorted(chromDict.keys()[:])
	if statsFlag:
		fig, axes = plt.subplots(nrows=len(sortedChroms)+1)
		fig.subplots_adjust(top=0.95, bottom=0.05, left=0.07, right=0.97)
		axes[0].set_title("RAT Heatmap")
		# Calculate the maximum differential in the heatmap
		numberTimes = len(nameList)
		distMax=numberTimes-1
	for chrom, ax in zip(sortedChroms, axes):
		chromMA = A[chrom]
		chromMB = B[chrom]
		dists = map(dist, chromMA, chromMB)
		#disps = map(abs, dists)
		if statsFlag:
			if len(chromMA) < 600:
				Y = np.array([np.nanmean(i) for i in grouper(dists, int(np.ceil(len(dists)/30.0)), fillvalue=0.0)])
			else:
				Y = np.array([np.nanmean(i) for i in grouper(dists, int(np.ceil(len(chromMA)/600.0)), fillvalue=0.0)])
			Y = np.vstack((Y,Y))
			ax.imshow(Y, aspect='auto', cmap=plt.get_cmap("RdBu"), interpolation='nearest', vmin=-distMax, vmax=distMax)
			pos = list(ax.get_position().bounds)
			x_text = pos[0]-0.01
			y_text = pos[1] + pos[3]/2.0
			fig.text(x_text, y_text, chrom, va='center', ha='right', fontsize=10)
			ax.set_axis_off()
		for index in np.where(np.absolute(dists) >= minD)[0]:
			s = index*tileSize
			e = s+tileSize
			strA = toSTR(chromMA[index])
			strB = toSTR(chromMB[index])
			imA = indexMean(chromMA[index])
			imB = indexMean(chromMB[index])
			yield("%s\t%i\t%i\t%.2f\t%s\t%s\t%.1f\t%.1f"%(chrom, s, e, dists[index], strA, strB, imA, imB))
	if statsFlag:
		fig.text(0.03, 0.5, "Chromosome", va='center', ha='center', rotation='vertical')
		cb1 = matplotlib.colorbar.ColorbarBase(axes[-1], cmap=plt.get_cmap("RdBu"), norm=matplotlib.colors.Normalize(vmin=-distMax, vmax=distMax), orientation='horizontal')
		plt.savefig("RAT Plot."+figExt)

def indexMean(BA):
	if not np.any(BA):
		return -1
	return np.mean(np.where(BA))

def grouper(iterable, n, fillvalue=None):
	"Collect data into fixed-length chunks or blocks"
	# grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
	args = [iter(iterable)] * n
	return izip_longest(*args, fillvalue=fillvalue)

def processGenome(chromDict, tileSize, gff, statsFlag, figExt):
	genome = makeGenomeStruct(chromDict, tileSize)
	updateGenomeStruct(genome, gff, tileSize, chromDict, statsFlag, figExt)
	return genome

def updateGenomeStruct(genome, gff, tileSize, chromDict, statsFlag, figExt):
	def process(location, name, tileSize):
		chrom, start, end = location
		binArray = toBA(name)
		sI = int(np.ceil(start/tileSize))
		eI = int(np.ceil(end/float(tileSize)))
		return (chrom, binArray, sI, eI)
	if statsFlag:
		segments = {chrom:[] for chrom in chromDict.keys()}
		for location, name in fileReader(gff):
			chrom, binArray, sI, eI = process(location, name, tileSize)
			if chrom in chromDict:
				genome[chrom][sI:eI] = binArray
				arrayStr =  ''.join(map(lambda x: str(int(x)), binArray))
				segments[chrom].append((arrayStr, eI-sI))
		title = os.path.splitext(os.path.split(gff)[1])[0]
		plotSize(segments, title, tileSize, figExt)
		plotComp(segments, title, tileSize, chromDict, figExt)
	else:
		for location, name in fileReader(gff):
			chrom, binArray, sI, eI = process(location, name, tileSize)
			if chrom in chromDict:
				genome[chrom][sI:eI] = binArray

def toBA(name):
	'''
	>>> toBA('ESMS')
	array([ True,  True, False], dtype=bool)
	>>> toBA('MSLS')
	array([False,  True,  True], dtype=bool)
	'''
	return np.array([N in name for N in nameList], dtype=np.bool)

def toSTR(BA):
	'''
	>>> toBA('ESMS')
	array([ True,  True, False], dtype=bool)
	>>> toBA('MSLS')
	array([False,  True,  True], dtype=bool)
	'''
	if not np.any(BA):
		return 'NA'
	return ''.join(compress(nameList, BA))

def plotVars(names):
	global nameList
	nameList = names.split(',')
	global colors
	if names == 'ES,MS,LS':
		#global times
		#times = ('ES','ESMS','MS','MSLS','LS','ESLS','ESMSLS')
		#myColors = ("#2250F1","#28C5CC","#1A8A12","#FFFD33","#FB0018","#EA3CF2","#FAB427")
		colors = ["#FB0018","#1A8A12","#FFFD33","#2250F1","#EA3CF2","#28C5CC","#FAB427"]
	else:
		colors = cm.rainbow(np.linspace(0,1,2**len(nameList)-1))

def plotComp(segments, title, tileSize, chromDict, figExt):
	plt.figure()
	yIndex = 0.1
	yHeight = 0.8
	sortedChroms = sorted(chromDict.keys())
	labels, inds, cinds = makeLabels()
	OT = open("composition_%s.tab"%(title), 'w')
	OT.write("Chr\t"+'\t \t'.join(labels)+'\t \tChr Length\n')
	for chrom in sortedChroms:
		otStr = '%s\t'%(chrom)
		chromSize = chromDict[chrom]
		X = np.zeros(2**len(nameList)-1)
		for arrayStr, size in segments[chrom]:
			sortedInd = inds[int(arrayStr,2)-1]
			X[sortedInd] += size*tileSize
		percents = list(np.round(X/float(chromSize),3))
		sP = map(lambda x: str(x*100)+'%', percents)
		otStr += '\t'.join([str(val) for tup in zip(X,sP) for val in tup])+'\t'+str(chromSize)+'\n'
		OT.write(otStr)
		xranges = zip(np.cumsum([0]+percents[:-1]), percents)
		plt.broken_barh(xranges, (yIndex, yHeight), lw=0, color=[colors[i] for i in cinds])
		yIndex += 1
	OT.close()
	plt.xlim((0,1))
	plt.yticks(np.arange(0.5, len(sortedChroms)), sortedChroms)
	plt.ylabel("Chromosome")
	plt.xlabel("Fraction of Chromosome")
	plt.title(title+" Chromosome Composition")
	patches = [mpatches.Patch(color=colors[cinds[i]], label=labels[i]) for i in xrange(len(labels))]
	plt.figlegend(patches, labels, loc='center right', ncol=1, frameon=False)
	plt.tight_layout(rect=[0,0,0.81,1.0])
	plt.savefig("composition_%s.%s"%(title, figExt))
	plt.close()

def makeLabels():
	'''
	>>> makeLabels()
	['ES', 'ESMS', 'ESLS', 'ESMSLS', 'MS', 'MSLS', 'LS']
	'''
	labels = []
	numNames = 2**len(nameList)-1
	for i in range(numNames):
		binRep = map(int, np.binary_repr(i+1,len(nameList)))
		name = ''
		for binI in range(len(binRep)):
			if binRep[binI]: name += nameList[binI]
		boolA = np.array(binRep, dtype=np.bool)
		val = np.mean(np.where(boolA))
		labels.append((name,val))
	sortedLabels = sorted(labels, key=itemgetter(1,0))
	inds = [sortedLabels.index(x) for x in labels]
	cinds = [labels.index(x) for x in sortedLabels]
	return (map(lambda x: x[0], sortedLabels), inds, cinds)

def plotSize(segments, title, tileSize, figExt):
	X = [[] for i in xrange(2**len(nameList)-1)]
	labels, inds, cinds = makeLabels()
	for chromList in segments.itervalues():
		for arrayStr, size in chromList:
			base10 = int(arrayStr,2)
			sortedInd = inds[base10-1]
			X[sortedInd].append(size*tileSize)
	print "%s Size Distribution"%(title)
	print "%-6s %10s %10s %10s %10s %10s %10s"%("","min","1st-Q","median","3rd-Q","max",'count')
	for segment, xIndex in zip(labels, range(len(labels))):
		try:
			fiveSum = fivenum(X[xIndex]) # (min, 1st-Q, median, 3rd-Q, max)
			args = (segment,)+fiveSum+(len(X[xIndex]),)
		except:
			args = (segment,)+(0,0,0,0,0,0)
		print "%-6s %10.1f %10.1f %10.1f %10.1f %10.1f %10i"%args
	plt.figure()
	plt.boxplot(X, labels=labels, showfliers=False)
	plt.ylabel("Segment Size (bp)")
	plt.xlabel("Time")
	plt.title(title+" Size Distribution")
	plt.savefig("size_dist_%s.%s"%(title, figExt))
	plt.close()
		
def fileReader(a):
	if not os.path.splitext(a)[1] == '.gff3':
		sys.exit("%s is not a gff3 file"%(a))
	for line in open(a,'r'):
		if line[0] != '#':
			yield(lineParser(line)) #((chrom, start, end), name)

# pre-compiled RE for finding the name in the GFF
nameRE = re.compile(r'Name=([^;]+);')

def lineParser(line):
	tmp = line.split('\t')
	try:
		location = (tmp[0], int(tmp[3])-1, int(tmp[4])) # (chrom, start, end)
	except:
		print "Couldn't parse:", tmp
		sys.exit()
	name = nameRE.search(tmp[8]).group(1) # name
	return location, name

def makeGenomeStruct(chromDict, tileSize):
	genome = {}
	for chrom, chromLen in chromDict.iteritems():
		numBins = int(np.ceil(chromLen/tileSize))
		genome[chrom] = np.zeros((numBins, 3), dtype=np.bool)
	return genome

def readFAI(inFile):
	'''
	Returns length of each chromosome in a dictionary.
	'''
	lineList = map(lambda x: x.rstrip('\n').split('\t'), open(inFile,'r').readlines())
	chromSizePairs = map(lambda x: (x[0],int(x[1])), lineList)
	return dict(chromSizePairs)

def dist(a,b):
	'''
	Calculates the shift of the index average between two binary arrays. If
	either of the arrays has no replication (all zero), the distance is
	returned as zero.

	   EML  index-mean	returned distance
	A: 110	0.5		1.5 - 0.5 = 1
	B: 011	1.5

	A: 001	2		1 - 2 = -1
	B: 010	1

	If the mean doesn't change, a distance of zero is returned instead. This will 
	only happen when moving between M, EL, and EML. This should be ok due to the
	low occurance of EL and EML in the results.

	   EML  index-mean	returned distance
	A: 101	1		0
	B: 010	1

	
	>>> dist(toBA('ESMS'),toBA('MSLS'))
	1.0
	>>> dist(toBA('LS'),toBA('MS'))
	-1.0
	>>> dist(toBA('ESLS'),toBA('MS'))
	0.0
	'''
	if np.sum(a) == 0 or np.sum(b) == 0:
		return 0.0
	indexMeanA = np.mean(np.where(a))
	indexMeanB = np.mean(np.where(b))
	return indexMeanB-indexMeanA

def fivenum(v):
	'''
	Returns Tukey's five number summary

	(min, 1st-Q, median, 3rd-Q, max)

	for the input vector, a list or array of numbers based on 1.5 times
	the interquartile distance
	'''
	try:
		np.sum(v)
	except TypeError:
		print('Error: you must provide a list or array of only numbers')
	naV = np.array(v)
	notNAN = np.logical_not(np.isnan(naV))
	q1 = np.percentile(naV[notNAN],25)
	q3 = np.percentile(naV[notNAN],75)
	#iqd = q3-q1
	md = np.median(naV[notNAN])
	#whisker = 1.5*iqd
	return np.min(naV[notNAN]), q1, md, q3, np.max(naV[notNAN])

if __name__ == "__main__":
	main()
else:
	nameList = ['ES','MS','LS']
