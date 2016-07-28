#!/usr/bin/env python

import numpy as np
from collections import Counter
import argparse, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from itertools import izip_longest
from operator import itemgetter

def main():
	parser = argparse.ArgumentParser(description="Investigates the flanking regiosn around EL and EML.")
	parser.add_argument("GFF3",metavar="GFF3", help="Segmentation Profile (mitotic)", nargs=1)
	parser.add_argument("-T",metavar="STR", help="Times (Default: %(default)s)", default="ES,MS,LS", type=str)
	args = parser.parse_args()
	plotVars(args.T)
	firstLast, allTimes = processGenome(args.GFF3[0])
	plt.figure(figsize=(4,8))
	ax = plt.subplot(211)
	ax.xaxis.set_visible(False)
	firstLastFreqs = countFreqs(firstLast)
	writeFreqs(firstLastFreqs, 'EL_links.tab')
	plotFreqs(firstLastFreqs, 'Frequency of EL Links')
	ax = plt.subplot(212)
	ax.xaxis.set_visible(False)
	allFreqs = countFreqs(allTimes)
	writeFreqs(allFreqs, 'EML_links.tab')
	plotFreqs(allFreqs, 'Frequency of EML Links')
	plt.tight_layout()
	plt.savefig('relationships.png')

def plotFreqs(counter, name):
	labels, inds, cinds = makeLabels()
	colDict = {labels[i]:colors[cinds[i]] for i in range(len(labels))}
	for k,v in counter.iteritems():
		plt.plot([0,1], [v,v], 'k-')
		plt.scatter([0],[v], s=60, c=colDict[k[0]])
		plt.scatter([1],[v], s=60, c=colDict[k[2]])
		plt.text(1.2, v, '-'.join(k).translate(None,'S'), fontsize=9)
	plt.xlim((-0.25,1.8))
	plt.title(name+' (total=%i)'%(sum(counter.values())), y=1.05)

def countFreqs(inList):
	c = Counter()
	c.update(inList)
	return c

def writeFreqs(counter, outFile):
	OF = open(outFile,'w')
	for k,i in counter.iteritems():
		OF.write('%s\t%i\n'%('\t'.join(k),i))
	OF.close()

def processGenome(gff):
	flName = nameList[0]+nameList[-1]
	allName = ''.join(nameList)
	firstLast = []
	allTimes = []
	for f, s, t in fileReader(gff):
		if s[1] == flName:
			if checkDist(f,s,t):
				firstLast.append((f[1], s[1], t[1]))
		elif s[1] == allName:
			if checkDist(f,s,t):
				allTimes.append((f[1], s[1], t[1]))
	return (firstLast, allTimes)

def checkDist(f,s,t):
	ntf = f[0][2] == s[0][1]
	ntt = s[0][2] == t[0][1]
	return ntf and ntt

def toBA(name):
	'''
	>>> toBA('ESMS')
	array([ True,  True, False], dtype=bool)
	>>> toBA('MSLS')
	array([False,  True,  True], dtype=bool)
	'''
	return np.array([N in name for N in nameList], dtype=np.bool)

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
		fiveSum = fivenum(X[xIndex]) # (min, 1st-Q, median, 3rd-Q, max)
		args = (segment,)+fiveSum+(len(X[xIndex]),)
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
	IF = open(a,'r')
	# Skip header
	tmpLine = IF.readline()
	while tmpLine[0] == '#': tmpLine = IF.readline()
	sLine = lineParser(tmpLine) #((chrom, start, end), name)
	tLine = lineParser(IF.readline())
	for line in IF:
		if line[0] != '#':
			fLine = sLine
			sLine = tLine
			tLine = lineParser(line)
			yield((fLine, sLine, tLine))

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

if __name__ == "__main__":
	main()
else:
	nameList = ['ES','MS','LS']
