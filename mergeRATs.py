#!/usr/bin/env python

import numpy as np
import scipy.signal as ss
import argparse, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from itertools import groupby, izip, imap
from StringIO import StringIO
from collections import defaultdict
from multiprocessing import Pool

def main():
	parser = argparse.ArgumentParser(description="Merge RATs from RATrap output")
	parser.add_argument("-z",metavar="INT", help="Merge regions separated by [%(default)s] bp of zero diff", default=0, type=int)
	parser.add_argument("-d",metavar="FLOAT", help="Differences in the same +/- direction that are smaller than '-d' will be merged [%(default)s]", default=0.0, type=float)
	parser.add_argument("--same", action="store_true", help="Only merge regions of the same time classification")
	parser.add_argument("--greedy", action="store_true", help="Merge small regions into larger ones")
	parser.add_argument("bg", metavar="BEDG", help="RATrap.py output bedgraph", nargs=1)
	args = parser.parse_args()
	# Cannot use both same and greedy
	if not (args.same or args.greedy):
		syse.xit("Please choose either same or greedy")
	if args.same and args.greedy:
		syse.xit("Please choose either same or greedy, not both")
	# Open Bedgraph
	IF = open(args.bg[0], 'r')
	if args.same:
		mergeSame(IF, args.z, args.d)
	elif args.greedy:
		mergeGreedy(IF, args.z, args.d)
	else:
		sys.exit("Not sure how I got here")

def mergeGreedy(IF, regionDist, ratDist):
	r'''
	>>> IF = StringIO('#Chromosome start end distance A B index-m_A index-m_B\n1\t51000\t54000\t-0.50\tMS\tESMS\t1.0\t0.5\n1\t54000\t57000\t-0.50\tMS\tESMS\t1.0\t0.5\n1\t111000\t114000\t-0.50\tESMS\tES\t0.5\t0.0\n2\t114000\t117000\t-0.50\tESMS\tES\t0.5\t0.0')
	>>> mergeSame(IF, 0, 0.5) # doctest: +NORMALIZE_WHITESPACE
	#Chromosome start end distance A B index-m_A index-m_B
	1	51000	57000	-0.5	MS	ESMS	1.0	0.5
	1	111000	114000	-0.5	ESMS	ES	0.5	0.0
	2	114000	117000	-0.5	ESMS	ES	0.5	0.0
	'''
	recA = ''
	records = imap(lambda x: tuple(x[:4]), readRec(IF))
	outRec = defaultdict(list)
	# Merge neighbor same
	for recB in records:
		if not recA:
			recA = recB
			continue
		if withinBound(recA, recB, regionDist, 0):
			recA = mergeRec(recA, recB)
			continue
		outRec[recA[0]].append(recA)
		recA = recB
	outRec[recA[0]].append(recA)
	del records
	# Merge starting from smallest
	p = Pool(initializer=initWorker((regionDist, ratDist)))
	chroms = sorted(outRec.keys())
	for chromRec in p.imap(mergeWorker, [outRec[c] for c in chroms]):
		for record in chromRec:
			printRec(record)

def initWorker(args):
	global regionDist, ratDist
	regionDist, ratDist = args
	
def mergeWorker(chromRec):
	modified = True
	while modified:
		for index in getSizeIndex(chromRec):
			modified = tryMerge(chromRec, index, regionDist, ratDist)
			if modified: break
	return chromRec
	
	

def tryMerge(records, index, regionDist, ratDist):
	'''
	>>> A = [(1, 0, 3, 2), (1, 3, 5, 1.5), (1, 5, 8, 1)]
	>>> tryMerge(A, 1, 0, 0.5)
	True
	>>> A
	[(1, 0, 4, 2), (1, 4, 8, 1)]
	>>> A = [(1, 0, 3, 2), (1, 3, 5, 1.5), (1, 5, 8, 1)]
	>>> tryMerge(A, 0, 0, 0.5)
	True
	>>> A
	[(1, 0, 5, 1.5), (1, 5, 8, 1)]
	>>> A = [(1, 0, 3, 2), (1, 3, 5, 1.5), (1, 5, 8, 1)]
	>>> tryMerge(A, 2, 0, 0.5)
	True
	>>> A
	[(1, 0, 3, 2), (1, 3, 8, 1.5)]
	>>> A = [(1, 0, 2, 2), (1, 3, 5, 1.5), (1, 5, 8, 1)]
	>>> tryMerge(A, 1, 0, 0.5)
	True
	>>> A
	[(1, 0, 2, 2), (1, 3, 8, 1)]
	>>> A = [(1, 0, 2, 2), (1, 3, 5, 1.5), (1, 5, 8, 1)]
	>>> tryMerge(A, 1, 1, 0.5)
	True
	>>> A
	[(1, 0, 3, 2), (1, 3, 8, 1)]
	'''
	target = records[index]
	# Check for edges
	if index > 0:
		leftRec = records[index-1]
		left = withinBound(target, leftRec, regionDist, ratDist)
	else:
		left = False
	if index < len(records)-1:
		rightRec = records[index+1]
		right = withinBound(target, rightRec, regionDist, ratDist)
	else:
		right = False
	if left and right:
		# split and merge both
		mergeThree(records, index-1, index, index+1)
	elif left:
		mergeTwo(records, index-1, index)
	elif right:
		mergeTwo(records, index+1, index)
	else:
		return False
	return True

def mergeTwo(A, a, b):
	'''
	Keep the value at index a, but merge in region b

	>>> A = [(1, 0, 3, 2), (1, 3, 5, 1.5), (1, 5, 8, 1)]
	>>> mergeTwo(A, 0, 1)
	>>> A
	[(1, 0, 5, 2), (1, 5, 8, 1)]
	>>> mergeTwo(A, 1, 0)
	>>> A
	[(1, 0, 8, 1)]
	'''
	tmpA = list(A[a])
	tmpB = list(A[b])
	assert tmpA[0] == tmpB[0]
	if tmpA[1] < tmpB[1]:
		tmpA[2] = tmpB[2]
	else:
		tmpA[1] = tmpB[1]
	A[a] = tuple(tmpA)
	del A[b]

def mergeThree(A, a, b, c):
	'''
	Split and merge b in a and c
	a gets floor len(b)/2
	c gets ceil len(b)/2

	>>> A = [(1, 0, 3, 2), (1, 3, 5, 1.5), (1, 5, 8, 1)]
	>>> mergeThree(A, 0, 1, 2)
	>>> A
	[(1, 0, 4, 2), (1, 4, 8, 1)]
	'''
	tmpA = list(A[a])
	tmpB = list(A[b])
	tmpC = list(A[c])
	assert tmpA[0] == tmpB[0] and tmpB[0] == tmpC[0]
	assert tmpA[1] < tmpB[1] and tmpB[1] < tmpC[1]
	lenB = tmpC[1]-tmpA[2]
	forA = int(lenB/2)
	forC = lenB-forA
	# Modify records
	tmpA[2] += forA
	tmpC[1] -= forC
	A[a] = tuple(tmpA)
	A[c] = tuple(tmpC)
	del A[b]

def getSizeIndex(records):
	'''
	Returns the index of records, sorted by record size

	>>> getSizeIndex([(1,0,3), (1,3,5), (1,5,9)])
	[1, 0, 2]
	'''
	return sorted(range(len(records)),key=lambda x:records[x][2]-records[x][1])

def toRegionIndex(A):
	'''
	Returns a region (value, size) index

	>>> toRegionIndex([2,2,2,1.5,1.5,1,1,1])
	[(2, 3, 0, 3), (1.5, 2, 3, 5), (1, 3, 5, 8)]
	'''
	index = 0
	OUT = []
	for value, t in groupby(A):
		size = len(tuple(t))
		OUT.append((value, size, index, index+size))
		index += size
	return OUT

def getDistofD(D):
	'''
	Returns a numpy array of distances
	'''
	return np.array([ R[3] for R in D ], dtype=np.dtype('f2'))

def sameClass(recA, recB):
	'''
	>>> sameClass([0,0,0,1,'ES','MS'],[0,0,0,1,'ES','MS'])
	True
	>>> sameClass([0,0,0,1,'ES','MS'],[0,0,0,1,'ESMS','MS'])
	False
	'''
	return (recA[4] == recB[4]) and (recA[5] == recB[5])
def withinBound(recA, recB, regionDist, ratDist):
	'''
	>>> withinBound([0,0,0,1],[0,0,0,1], 0, 0)
	True
	>>> withinBound([0,0,1,1],[0,1,2,1], 0, 0)
	True
	>>> withinBound([0,0,1,1],[0,1,2,1.5], 0, 0.5)
	True
	>>> withinBound([0,0,1,1],[0,1,2,2], 0, 0.5)
	False
	>>> withinBound([0,0,1,1],[0,2,3,1], 0, 0)
	False
	>>> withinBound([1,0,0,1],[0,0,0,1], 0, 0)
	False
	>>> withinBound([0,0,0,1],[0,0,0,-1], 0, 0)
	False
	'''
	sameChrom = recA[0] == recB[0]
	oDist = abs(recA[2]-recB[1])
	ooDist = abs(recB[2]-recA[1])
	ltRegionDist = oDist <= regionDist or ooDist <= regionDist
	ltRatDist = abs(recA[3] - recB[3]) <= ratDist
	return sameChrom and ltRegionDist and ltRatDist
def mergeRec(recA, recB):
	#Chromosome	start	end	distance	A	B	index-m_A	index-m_B
	#1	44000	45000	-0.50	MS	ESMS	1.0	0.5
	#1	45000	46000	-0.50	MS	ESMS	1.0	0.5
	assert recA[0] == recB[0]
	assert recA[2] <= recB[1]
	tmp = list(recA)
	tmp[2] = recB[2]
	return tuple(tmp)
def mergeSame(IF, regionDist, ratDist):
	r'''
	>>> IF = StringIO('#Chromosome start end distance A B index-m_A index-m_B\n1\t51000\t54000\t-0.50\tMS\tESMS\t1.0\t0.5\n1\t54000\t57000\t-0.50\tMS\tESMS\t1.0\t0.5\n1\t111000\t114000\t-0.50\tESMS\tES\t0.5\t0.0\n2\t114000\t117000\t-0.50\tESMS\tES\t0.5\t0.0')
	>>> mergeSame(IF, 0, 0.5) # doctest: +NORMALIZE_WHITESPACE
	#Chromosome start end distance A B index-m_A index-m_B
	1	51000	57000	-0.5	MS	ESMS	1.0	0.5
	1	111000	114000	-0.5	ESMS	ES	0.5	0.0
	2	114000	117000	-0.5	ESMS	ES	0.5	0.0
	'''
	recA = ''
	for recB in readRec(IF):
		if not recA:
			recA = recB
			continue
		if withinBound(recA, recB, regionDist, ratDist):
			if sameClass(recA, recB):
				recA = mergeRec(recA, recB)
				continue
		printRec(recA)
		recA = recB
	printRec(recA)

def printRec(rec):
	if len(rec) == 8:
		print "%s\t%i\t%i\t%0.1f\t%s\t%s\t%0.1f\t%0.1f"%tuple(rec)
	elif len(rec) == 4:
		print "%s\t%i\t%i\t%0.1f"%tuple(rec)
	else:
		sys.exit("Unhandled record length for printing")

def readRec(FP):
	r'''
	>>> TF = StringIO("#Chromosome start end distance A B index-m_A index-m_B\n1\t44000\t45000\t-0.50\tMS\tESMS\t1.0\t0.5\n1\t45000\t46000\t-0.50\tMS\tESMS\t1.0\t0.5")
	>>> for i in readRec(TF): print i
	#Chromosome start end distance A B index-m_A index-m_B
	('1', 44000, 45000, -0.5, 'MS', 'ESMS', 1.0, 0.5)
	('1', 45000, 46000, -0.5, 'MS', 'ESMS', 1.0, 0.5)
	'''
	for line in FP:
		if line[0] == "#":
			print line.rstrip('\n')
			continue
		tmp = line.rstrip('\n').split('\t')
		yield (tmp[0], int(tmp[1]), int(tmp[2]), float(tmp[3]), tmp[4], tmp[5], float(tmp[6]), float(tmp[7]))


if __name__ == "__main__":
	main()
