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

## Running the Pipeline

### Usage
`logFC_calculate.py [-h] -F FASTA [-L INT] [-S INT] FILE`

### Required Arguments

| Flag | Option | Description |
|:----:|:------:|-------------|
|-F|FASTA|The fasta file used for alignment|
|-L|INT|The level of smoothing to use [1,5] (Default: 2)|
|-S|INT|The size of each window in the bedgraphs (Default: 500)|
|  |TXT| A text file listing bams for input|

### Input TXT
Each line of the text file needs to contain a path to a bam file and a short name for the file separated by a **tab**. For example:
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

## Output
| File | Description |
|:----:|-------------|
|`*.bedgraph`|Bedgraph produced from corresponding bam|
|`*_logFC.bedgraph`|Bedgraph of signal after dividing control and performing log2 transform|
|`*_logFC_*.smooth.bedgraph`|Smoothed using specified level of Haar wavelet|
|`*_logFC_*.smooth.gff3`|GFF showing positive regions|
|`logFC_segmentation.gff3`| Segmentation GFF|
