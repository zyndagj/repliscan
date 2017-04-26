#!/usr/bin/python

"""
Setup script for bioitools
"""
try:
	from setuptools import setup, Extension
except:
	from distutils.core import setup, Extension

# Make wavelets
import subprocess
subprocess.call(['make'])

setup(name = "repliscan",
	version = "1.0.0",
	author = "Greg Zynda",
	author_email="gzynda@tacc.utexas.edu",
	license="GNU",
	description="Pipeline for calculating Repli-seq enrichment and classifying the time replication took place.",
	scripts = ["repliscan.py","RATrap.py","wavelets"])
