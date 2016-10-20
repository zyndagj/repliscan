# Repliscan - Logfold Reptiming Pipeline
Pipeline for calculating reptiming enrichment and classifying the time replication took place.

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
4. python

   From Stampede @ TACC
   ```
   $ module load python
   ```

## Methods

![Workflow DAG](dag.jpg)

## Running the Pipeline

### Usage
```
repliscan.py [-h] -r FASTA [-l INT] [-w INT] [-a STR] [-n STR] [-t STR]
             [-S STR] [-v Float] [-p Float] [-c STR] [-R STR] [--log]
             [-f] [--plot] FILE
```

### Input TXT - FILE
Each line of the text file needs to contain a short name describing the sample and then a list of bam files corresponding to that name, all separated by *tabs*.

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

### Arguments

| Flag | Option | Description - Bold denotes Default|
|:----:|:------:|:------------|
|-r, --ref|**FASTA**|Reference fasta|
|-l,--level|INT|Haar smoothing level \[1,2,**3**,4,5\]|
|-w, --window|INT|Analysis bin size in base pairs - **1000**|
|-a, --aggregate|STR|Replicate agregation method \(**sum**, median, mean, min, max\)|
|-n, --norm|STR|Normalization Method \[DESeq, **Coverage**\]|
|-t, --threshold|STR|Replication threshold method \[value|**auto**|percent\]|
|-S,--scope|STR|Replication scope \[**chromosome**|genome\]|
|-v, --value|Float|Explicit replication threshold value \[1.0\]|
|-p, --percent|Float|Replication percent cut \[2.0\]|
|-c, --classifier|STR|Segmentation classification method \[binary|**proportion**\]|
|-R, --remove|STR|Outlying data to remove \[none|gamma|**norm**|whiskers\]|
|--log| |Apply log transform to sequenceability ratio|
|-f, --force| |Force the re-generation of all files|
|--plot| |Plot Statistics|
|  |FILE| A text file listing bams for input ***Required***|

### Sequencability Method
- `default` - \(sample/control\)
- `--log` - log\(sample/control\) *Default*

### Normalization Methods
- `DESeq` - DESeq size normalization using geometric mean.
- `Coverage` - Transform each sample to 1X coverage across the genome. This replicates the RPGC normalization method from deepTools.

### Replication Methods
- `value` - A signal ratio above threshold `-v` is considered replicating.
- `auto` - Determine a signal ratio threshold unique to each chromosome based on change in coverage.
- `percent` - Signal ratio values for each chromosome above percentile `-p` are considered replicating.

### Segmentation Methods
- `binary` - Time classifications are combined on a binary basis.
- `proportion` - Time classifications are determined based on proportion.

### Handling Replicates
If you have replicates in your input file,
```
MS	MS_001.bam	MS_L1.bam	MS_L2.bam
```
you can specify the method by which they are aggregated with the `-a` paramter. This accepts the following methods:
  - sum (Default)                                    
  - median                                           
  - mean                                             
  - min                                              
  - max                                              

## Output
 - `[replicate].bed` - Bed files produced from corresponding bam
 - `[replicate].bedgraph` - Bedgraph produced from corresponding bed
 - `[replicate]_norm.bedgraph` - Normalized bedgraph produced from corresponding bed
 - `[time].bedgraph` - Aggregated signal from replicates
 - `[time]_norm.bedgraph` - Normalized aggregated signal from replicates
 - `[time]_(logFC|ratio).bedgraph` - Ratio signal when each time is normalized for sequenceability
 - `[time]_(logFC|ratio)_*.smooth.bedgraph` - Ratio signal smoothed using a specified level of Haar wavelet
 - `[time]_(logFC|ratio)_*.smooth.gff3` - GFF showing replicating regions
 - `(logFC|ratio)_segmentation.gff3` -  Segmentation GFF, which can be used as input for `RATrap.py`
