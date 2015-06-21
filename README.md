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
`python logfold_rep.py [-h] -F FASTA [-L INT] [-S INT] [-C STR] FILE`

### Required Arguments

| Flag | Option | Description |
|:----:|:------:|:------------|
|-F|FASTA|The fasta file used for alignment|
|-L|INT|The level of smoothing to use \[1,5\] \(Default: 2\)|
|-S|INT|The size of each window in the bedgraphs \(Default: 500\)|
|-C|STR|How to handle replicates \(Default: sum\)|
|--rep STR|Replicating Method \(threshold|auto|percent\) \(Default: threshold\)|
|--seg STR|Segmentation Method \(binary|proportion\) \(Default: binary\)|
|-T Float|Threshold Level (Default: 0.0)|
|-P Float|Percent Cut \(Default: 2.0\)|
|--plot|Plot Coverage|
|  |TXT| A text file listing bams for input|

### Input TXT
Each line of the text file needs to contain a short name describing the sample and then a list of bam files corresponding to that name, all separated by tabs.

For example:

```
EX	example_rep_1.bam	example_rep_2.bam
```
The first line of this file needs to be the control (G1). All subsequent lines need to be listed sequentially according to experimental time. An example file would be:
```
G1	G1_001.bam	G1_002.bam
ES	ES_001.bam
MS	MS_001.bam	MS_L1.bam	MS_L2.bam
LS	LS.bam
```

### Replication Methods
- threshold - A log(ratio) above this value is considered to be replicating.
- auto - Determine a log(ratio) threshold unique to each chromosome based on change in coverage.
- percent - Values greater than this percentile are considered to be replicating.

### Segmentation Methods
- binary - Time classifications are combined on a binary basis.
- proportion - Time classifications are determined based on proportion.

### Handling Replicates
  - sum (Default)                                    
  - median                                           
  - mean                                             
  - min                                              
  - max                                              

## Output
| File | Description |
|:----:|-------------|
|`*.bedgraph`|Bedgraph produced from corresponding bam (not normalized)|
|`*_logFC.bedgraph`|Bedgraph of signal after dividing control and performing log2 transform|
|`*_logFC_*.smooth.bedgraph`|Smoothed using specified level of Haar wavelet|
|`*_logFC_*.smooth.gff3`|GFF showing positive regions|
|`logFC_segmentation.gff3`| Segmentation GFF|
