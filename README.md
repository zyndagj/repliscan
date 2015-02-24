# Logfold Reptiming Pipeline
Pipeline for calculating reptiming enrichment using logfold and Haar wavelet smoothing.
## Dependencies
The following binaries need to exist on the user's PATH:

1. bedtools - https://github.com/arq5x/bedtools2
   
   From Stampede/Lonestar @ TACC
   ```
   $ module load bedtools
   ```
2. samtools - https://github.com/samtools/samtools
   
   From Stampede/Lonestar @ TACC
   ```
   $ module load samtools
   ```
3. wavelets - http://staff.washington.edu/dbp/WMTSA/NEPH/wavelets.html
   
   Installation
   ```
   $ tar -xzf wavelets.tgz
   $ cd wavelets/src
   $ make
   $ cd ../bin
   $ cp wavelets [to somewhere on bin path]
   ```
## Running `logfold_rep.py`

`usage: logFC_calculate.py [-h] -F FASTA [-L INT] [-S INT] FILE`

### Required Arguments

| Flag | Option | Description |
|:----:|--------|-------------|
|-F|Fasta File|The fasta file used for alignment|
|-L|Smoothing Level|The level of smoothing to use [1,5]|
|-S|Window Size|The size of each window in the bedgraphs|
|  |BAM List File| A text file listing bams for input|

### Bam List File
Each line of the bam list file needs to contain a path to a bam file and a short name for the file separated by a **tab**. For example:
```
example_sorted.bam	EX
```
The first line of this file needs to be the control (G1). All subsequent lines need to be listed sequentially according to experimental time. An example file would be:
```
G1_sorted.bam	G1
ES_sorted.bam	ES
MS_sorted.bam	MS
LS_sorted.bam	LS
```
