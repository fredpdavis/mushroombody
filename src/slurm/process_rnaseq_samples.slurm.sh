#!/bin/csh -v
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-core=1
#SBATCH --mem=20g
#SBATCH -o slurm_out/process_rnaseq_samples.%A.%a.out
#SBATCH -e slurm_out/process_rnaseq_samples.%A.%a.err
#SBATCH --time=10:00:00
#SBATCH --gres=lscratch:200
#SBATCH --array 1-14


set NUMCPU=4

set BASEDIR="/data/davisfp/projects/dubnau.mb"
set SAMPLEINFO_FN="$BASEDIR/meta/kenyoncell_ms_rnaseq_samples.txt"

set samplenums=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN fastqSampleNum` )
set samplenames=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN sampleName` )
set SAMPLENUM=$samplenums[$SLURM_ARRAY_TASK_ID]
set SAMPLENAME=$samplenames[$SLURM_ARRAY_TASK_ID]

## LOAD PROGRAMS
module load seqtk/1.2
module load STAR/2.5.3a
module load picard/1.119
module load java/1.8.0_92
module load R/3.3.2_gcc-4.9.1
module load deeptools/2.5.0.1

set KALLISTO_BIN="~/software/kallisto/kallisto_linux-v0.43.1/kallisto"
set PICARDSORT_BIN="java -jar $PICARD_JARPATH/SortSam.jar"
set PICARDINDEX_BIN="java -jar $PICARD_JARPATH/BuildBamIndex.jar"
set PICARDSTATS_BIN="java -jar $PICARD_JARPATH/CollectRnaSeqMetrics.jar"


set trimlength=5


## INPUT FILES
set FASTQ1_FN="$BASEDIR/data/fastq/S${SAMPLENUM}_R1_001.fastq.gz"
set FASTQ2_FN="$BASEDIR/data/fastq/S${SAMPLENUM}_R2_001.fastq.gz"

set KALLISTO_INDEX="$BASEDIR/data/kallisto_files.BDGP6.91/BDGP6.91.ERCC.INTACT.ncrna.kallisto_index"

set STAR_INDEX="$BASEDIR/data/star_files.BDGP6.91/BDGP6.91.ERCC.INTACT.star_index"

set GENOME_FASTAF="$BASEDIR/data/star_files.BDGP6.91/BDGP6.91.ERCC.INTACT.fa"
set STAR_GTF="$BASEDIR/data/star_files.BDGP6.91/BDGP6.91.ERCC.INTACT.gtf"
set REFFLAT_FN="$BASEDIR/data/picard_files.BDGP6.91/refFlat.noERCC.txt.gz"
set ERCC_REFFLAT_FN="$BASEDIR/data/picard_files.BDGP6.91/refFlat.ERCConly.txt.gz"


## OUTPUT DIRECTORIES
set TRIM_FASTQF_DIR="$BASEDIR/results/RNAseq/trim_reads/$SAMPLENUM"
set KALLISTO_OUTDIR="$BASEDIR/results/RNAseq/kallisto.BDGP6.91/$SAMPLENUM"
set STAR_OUTDIR="$BASEDIR/results/RNAseq/star_align.BDGP6.91/$SAMPLENUM"
set PICARDSTATS_OUTDIR="${BASEDIR}/results/RNAseq/picard_stats.star.BDGP6.91/$SAMPLENUM"

foreach T_OUTDIR ( $TRIM_FASTQF_DIR $KALLISTO_OUTDIR $STAR_OUTDIR $PICARDSTATS_OUTDIR )
   if (! -e $T_OUTDIR) then
      echo "mkdir -p $T_OUTDIR"
      mkdir -p $T_OUTDIR
   endif
end


### OUTPUT FILES
set TRIM_FASTQ1_FN="$TRIM_FASTQF_DIR/${SAMPLENUM}_R1_trim.fastq.gz"

set STAR_OUT_BAM_FN="${STAR_OUTDIR}/Aligned.out.bam"
set STAR_OUT_SORTED_BAM_FN="${STAR_OUTDIR}/Aligned.out.sorted.bam"
set STAR_OUT_SORTED_BAI_FN="${STAR_OUTDIR}/Aligned.out.sorted.bai"
set STAR_OUT_BW_FN="${STAR_OUTDIR}/${SAMPLENUM}.star_deeptool_rpkm.bw"
set STAR_NAMED_OUT_BW_FN="${STAR_OUTDIR}/${SAMPLENAME}.star_deeptool_rpkm.bw"

set PICARDSTATS_TXT_FN="${PICARDSTATS_OUTDIR}/${SAMPLENUM}.noERCC.picard_rnaseq_report.txt"
set PICARDSTATS_PDF_FN="${PICARDSTATS_OUTDIR}/${SAMPLENUM}.noERCC.picard_rnaseq_report.pdf"
set PICARDSTATS_ERCC_TXT_FN="${PICARDSTATS_OUTDIR}/${SAMPLENUM}.ERCConly.picard_rnaseq_report.txt"
set PICARDSTATS_ERCC_PDF_FN="${PICARDSTATS_OUTDIR}/${SAMPLENUM}.ERCConly.picard_rnaseq_report.pdf"



### PROGRAM OPTIONS
set KALLISTO_OPTIONS="quant -i $KALLISTO_INDEX -b 50 -o $KALLISTO_OUTDIR $TRIM_FASTQ1_FN $FASTQ2_FN"
set STAR_OPTIONS="--genomeDir $STAR_INDEX --readFilesIn $TRIM_FASTQ1_FN $FASTQ2_FN --outSAMtype BAM Unsorted --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --runThreadN $NUMCPU --outFileNamePrefix $STAR_OUTDIR/ --readFilesCommand zcat"


set PICARDSORT_OPTIONS="I=$STAR_OUT_BAM_FN O=$STAR_OUT_SORTED_BAM_FN SO=coordinate"
set PICARDINDEX_OPTIONS="I=$STAR_OUT_SORTED_BAM_FN O=$STAR_OUT_SORTED_BAI_FN"
set PICARDSTATS_OPTIONS1="REF_FLAT=$REFFLAT_FN STRAND_SPECIFICITY=NONE INPUT=$STAR_OUT_SORTED_BAM_FN CHART_OUTPUT=$PICARDSTATS_PDF_FN OUTPUT=$PICARDSTATS_TXT_FN"
set PICARDSTATS_OPTIONS2="REF_FLAT=$ERCC_REFFLAT_FN STRAND_SPECIFICITY=NONE INPUT=$STAR_OUT_SORTED_BAM_FN CHART_OUTPUT=$PICARDSTATS_ERCC_PDF_FN OUTPUT=$PICARDSTATS_ERCC_TXT_FN"

set BAMCOVERAGE_OPTIONS="-b $STAR_OUT_SORTED_BAM_FN -p $NUMCPU -o $STAR_OUT_BW_FN --normalizeUsingRPKM"





### START ACTUALLY RUNNING

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`

echo "# slurm run started on $curhost at $curtime"
echo "# Processing sample $SAMPLENUM" ;

set echo

set curtime=`date`
echo "STEP 1. READ TRIMMING ($curtime)"
seqtk trimfq -b $trimlength ${FASTQ1_FN} | gzip > $TRIM_FASTQ1_FN


set curtime=`date`
echo "STEP 2. KALLISTO PSEUDOALIGNMENT ($curtime)"
$KALLISTO_BIN $KALLISTO_OPTIONS


set curtime=`date`
echo "STEP 3. STAR ($curtime)"
STAR $STAR_OPTIONS

set curtime=`date`
echo "STEP 4. PICARD ($curtime)"
$PICARDSORT_BIN $PICARDSORT_OPTIONS
$PICARDINDEX_BIN $PICARDINDEX_OPTIONS

bamCoverage $BAMCOVERAGE_OPTIONS
ln -s $STAR_OUT_BW_FN $STAR_NAMED_OUT_BW_FN

$PICARDSTATS_BIN $PICARDSTATS_OPTIONS1
$PICARDSTATS_BIN $PICARDSTATS_OPTIONS2

set curtime=`date`
echo "# slurm run finished on $curhost at $curtime"
