# mushroombody

This package contains code to analyze RNA-seq data and generate figures and
tables presented in:

Nuclear transcriptomes of the seven neuronal cell types that constitute the
Drosophila mushroom bodies. 
Shih MF*, Davis FP*, Henry GL+, Dubnau J+

Contact fred.davis@nih.gov with any questions.

## Package contents

- src - code, organized by language, to turn RNA-seq read files into figures.
- metadata - text tables describing RNA-seq samples
- data - text data files used by code

## Requirements

### Data

The data used by this package comes from several sources. We include nearly all
data files expected by the R and slurm shell programs, along with README files
describing the contents. The only exceptions are large files (eg, genome
sequence, gene annotations, transcript sequences, RNA-seq alignment indices),
which we do not provide but describe in README files how to obtain or build.

| Data                    | Source                                                                                                                 |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| RNA-seq FASTQ files     | [GEO accession GSE119629](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119629)                                |
| genome sequence         | [ENSEMBL release 91; based on BDGP6 (dm6)](http://dec2017.archive.ensembl.org/Drosophila_melanogaster/Info/Index)      |
| gene structures         | [ENSEMBL release 91; based on FlyBase 2017_04](http://dec2017.archive.ensembl.org/Drosophila_melanogaster/Info/Index)  |


### Software

We used the following software on the [NIH Biowulf](https://hpc.nih.gov)
slurm-based linux cluster.

|  #  | Software               |  Source                                    |
| --- | ---------------------- | ------------------------------------------ |
|  1  | seqtk 1.2              | https://github.com/lh3/seqtk               |
|  2  | kallisto 0.43.1        | https://pachterlab.github.io/kallisto/     |
|  3  | STAR 2.5.3a            | https://github.com/alexdobin/STAR          |
|  4  | picard 1.119           | http://broadinstitute.github.io/picard/    |
|  5  | R v3.5.0               | https://www.r-project.org/                 |
|  6  | deeptools 2.5.0.1      | https://deeptools.readthedocs.io/          |

## Analysis overview

The analysis includes procerssing the RNA-seq reads and then generating figures
and tables.

## 1. RNA-seq read processing

This step is implemented in a slurm script that trims the reads (seqtk),
pseudo-aligns them to the transcriptome (kallisto) and aligns them to the
genome (STAR), evaluates quality metrics (picard), and generates bigwig
tracks for visualization (samtools, ucsc kent tools)

```sh
sbatch process_rnaseq_samples.slurm.sh
```

## 2. Generating figures and tables

An R script generates all the tables and figures in the manuscript.

```R
source("../../src/R/analyzeKenyonTapinSeq_ms.R")
dat <- main()
```
