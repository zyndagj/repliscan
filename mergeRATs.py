#!/usr/bin/env python

import numpy as np
import argparse, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm

def main():
	parser = argparse.ArgumentParser(description="Merge RATs from RATrap output")
	parser.add_argument("-t",metavar="FLOAT", help="Differences in the same +/- direction that are smaller than '-t' will be merged [%(default)s]", default=0.0, type=float)
	parser.add_argument("--same", action="store_true", help="Only merge regions of the same time classification")
	parser.add_argument("bg", metavar="BEDG", help="RATrap.py output bedgraph")
	args = parser.parse_args()
	IF = open(args.bg, 'r')
	print IF.readline().rstrip('\n')
	recA = ''
	for line in readRec(IF):
		if not recA:
			recA = line
			continue
		# Same chromosome and adjacent
		if recA[0] == line[0] and recA[2] == line[1]:
			# Same sign
			if sameSign(recA, line):
				# Mean within bound
				if withinBound(recA, line, args.t):
					# Same class
					if args.same:
						if sameClass(recA, line):
							recA = mergeRec(recA, line)
						else:
							printRec(recA)
							recA=line
					else:
						recA = mergeRec(recA, line)
				else:
					printRec(recA)
					recA=line
			else:
				printRec(recA)
				recA=line
		else:
			printRec(recA)
			recA = line
	printRec(recA)

def sameClass(recA, recB):
	'''
	>>> sameClass([0,0,0,1,'ES','MS'],[0,0,0,1,'ES','MS'])
	True
	>>> sameClass([0,0,0,1,'ES','MS'],[0,0,0,1,'ESMS','MS'])
	False
	'''
	return (recA[4] == recB[4]) and (recA[5] == recB[5])

def withinBound(recA, recB, bound):
	'''
	>>> withinBound([0,0,0,1],[0,0,0,1], 0)
	True
	>>> withinBound([0,0,0,1],[0,0,0,-1], 0)
	False
	'''
	sA = recA[3]
	sB = recB[3]
	return abs(sA - sB) <= bound

def sameSign(recA, recB):
	'''
	>>> sameSign([0,0,0,1],[0,0,0,1])
	True
	>>> sameSign([0,0,0,1],[0,0,0,-1])
	False
	'''
	sA = recA[3]
	sB = recB[3]
	return (sA < 0) == (sB < 0)
		
def printRec(rec):
	if len(rec) == 8:
		print "%s\t%i\t%i\t%0.1f\t%s\t%s\t%0.1f\t%0.1f"%tuple(rec)
	elif len(rec) == 4:
		print "%s\t%i\t%i\t%0.1f"%tuple(rec)
	else:
		sys.exit("Unhandled record length for printing")

def readRec(FP):
	for line in FP:
		tmp = line.rstrip('\n').split('\t')
		yield [tmp[0], int(tmp[1]), int(tmp[2]), float(tmp[3]), tmp[4], tmp[5], float(tmp[6]), float(tmp[7])]

def mergeRec(recA, recB):
	#Chromosome	start	end	distance	A	B	index-m_A	index-m_B
	#1	44000	45000	-0.50	MS	ESMS	1.0	0.5
	#1	45000	46000	-0.50	MS	ESMS	1.0	0.5
	assert recA[0] == recB[0]
	assert recA[2] == recB[1]
	if recA[4] == recB[4] and recA[5] == recB[5]:
		tmp = recA
		tmp[2] = recB[2]
	else:
		tmp = recA[:4]
		# Take the index mean of the second record to merge thigns like
		# 0, 0.5, 1
		tmp[3] = recB[3]
	return tmp

if __name__ == "__main__":
	main()
