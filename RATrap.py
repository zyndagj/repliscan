#==> 07_seg.gff3 <==
###gff-version 3
##track name="Segmentation 1000bp" gffTags=on
#1	.	gene	401001	422000	.	.	.	ID=gene1;Name=LS;color=#FB0018;
#1	.	gene	809001	815000	.	.	.	ID=gene2;Name=LS;color=#FB0018;
#
#==> 09_seg.gff3 <==
###gff-version 3
##track name="Segmentation 1000bp" gffTags=on
#1	.	gene	236001	248000	.	.	.	ID=gene1;Name=LS;color=#FB0018;
#1	.	gene	329001	366000	.	.	.	ID=gene2;Name=LS;color=#FB0018;
#
import numpy as np
import argparse, os, re
from itertools import izip

# Didn't know how to handle EML and EL, so I'm treating distance as hamming to
# the inclusion
times = ('ES','ESLS','ESMS','ESMSLS','LS','MS','MSLS')
boolTimes = {'E':np.array((1,0,0),dtype=np.bool),'M':np.array((0,1,0),dtype=np.bool),'L':np.array((0,0,1),dtype=np.bool)}
timeRE = re.compile(r'Name=(([EML]S){1,3})')

def main():
	parser = argparse.ArgumentParser(description="Finds the timing differences between two segmentation profiles.")
	parser.add_argument("-d",metavar="INT", help="Minimum distance to be RAT (Default: %(default)s)", default=2, type=int)
	parser.add_argument("-S",metavar="INT", help="Tile Size (Default: %(default)s)", default=1000, type=int)
	parser.add_argument("-A",metavar="GFF3", help="Segmentation Profile", required=True)
	parser.add_argument("-B",metavar="GFF3", help="Segmentation Profile", required=True)
	parser.add_argument("-F",metavar="FASTA", help="Reference", required=True)
	parser.add_argument("-O",metavar="BEDG", help="Output to bedgraph file")
	args = parser.parse_args()
	if os.path.splitext(args.F)[1] in ['.fasta','.fa']:
		fai = args.F+'.fai'
	else:
		sys.exit("Please specify a fasta file\n")
	chromDict = readFAI(fai)
	genomeA = processGenome(chromDict, args.S, args.A)
	genomeB = processGenome(chromDict, args.S, args.B)
	if args.O:
		OF = open(args.O,'w')
		for record in compareGenomes(genomeA, genomeB, chromDict, args.d, args.S):
			OF.write(record+'\n')
		OF.close()
	else:
		for record in compareGenomes(genomeA, genomeB, chromDict, args.d, args.S):
			print record

def compareGenomes(A, B, chromDict, minD, tileSize):
	sortedChroms = sorted(chromDict.keys()[:])
	for chrom in sortedChroms:
		chromMA = A[chrom]
		chromMB = B[chrom]
		oldRec = False
		for index in xrange(chromMA.shape[0]):
			if dist(chromMA[index], chromMB[index]) >= minD:
				if oldRec:
					if index == oldRec[1]:
						oldRec[1] = index+1
					else:
						s = oldRec[0]*tileSize
						e = min(chromDict[chrom], oldRec[1]*tileSize)
						yield("%s\t%i\t%i"%(chrom, s,e))
						oldRec = [index, index+1]
				else:
					oldRec = [index, index+1]
		if oldRec:
			s = oldRec[0]*tileSize
			e = min(chromDict[chrom], oldRec[1]*tileSize)
			yield("%s\t%i\t%i"%(chrom, s,e))

def processGenome(chromDict, tileSize, gff):
	genome = makeGenomeStruct(chromDict, tileSize)
	updateGenomeStruct(genome, gff, tileSize)
	return genome

def updateGenomeStruct(genome, gff, tileSize):
	for location, color in fileReader(gff):
		chrom = location[0]
		binArray = toBA(color)
		sI = np.ceil(location[1]/tileSize)
		eI = np.ceil(location[2]/float(tileSize))
		genome[chrom][sI:eI] = binArray
		
def fileReader(a):
	if not os.path.splitext(a)[1] == '.gff3':
		sys.exit("%s is not a gff3 file"%(a))
	IF = open(a,'r')
	for line in IF:
		if line[0] != '#':
			yield(lineParser(line)) #((chrom, start, end), color)
	IF.close()

def lineParser(line):
	tmp = line.split('\t')
	location = (tmp[0], int(tmp[3])-1, int(tmp[4]))
	color = timeRE.search(tmp[8]).group(1)
	return location, color

def makeGenomeStruct(chromDict, tileSize):
	genome = {}
	for chrom, chromLen in chromDict.iteritems():
		genome[chrom] = makeChromArray(chromLen, tileSize)
	return genome

def makeChromArray(chromLen, tileSize):
	return np.zeros((np.ceil(chromLen/tileSize), 3), dtype=np.bool)

def readFAI(inFile):
	'''
	Returns length
	'''
	#Pt     140384  50      60      61
	chromDict = {}
	for line in open(inFile,'r'):
		tmp = line.split('\t')
		chromDict[tmp[0]] = int(tmp[1])
	return chromDict

def dist(a,b):
	'''
	>>> dist(toBA('ESMS'),toBA('MSLS'))
	2
	>>> dist(toBA('ES'),toBA('MS'))
	2
	>>> dist(toBA('ES'),toBA('ESMS'))
	1
	>>> dist(toBA('ES'),np.array([0,0,0]))
	0
	'''
	if np.sum(a) == 0 or np.sum(b) == 0:
		return 0
	return np.sum(np.logical_xor(a, b))

def toBA(a):
	'''
	>>> toBA('ESMS')
	array([ True,  True, False], dtype=bool)
	>>> toBA('MSLS')
	array([False,  True,  True], dtype=bool)
	'''
	return np.logical_or.reduce(map(lambda x: boolTimes[x], a.split('S')[:-1]))

if __name__ == "__main__":
	main()
