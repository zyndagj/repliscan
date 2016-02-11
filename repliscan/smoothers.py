#!/usr/bin/env python

import os, sys
import numpy as np
import pywt

def wavelets(level, fList, chromDict, use):
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
