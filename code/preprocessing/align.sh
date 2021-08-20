#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=100G
#$ -l mem_free=6G
#$ -l h_vmem=6G
#$ -l h_rt=96:00:00
#$ -pe local 8
#$ -t 1-41

module load bowtie2
module load samtools
CORES=8

umask g+w

# Inputs
fqdir="../fastq"
outdir="../bam"
scratchdir="tmp"
## align with alternates
refgenome="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/genome/hg19/hg19"
outdirfastp=$outdir/fastp

# Main
mkdir -p $outdir $outdir/dup_metrics $outdir/fastp $beddir $scratchdir
samplepath=$(find $fqdir -maxdepth 1 -name "*.fastq.gz" | \
    sed 's/.\{12\}$//' | \
    sort -u | \
    head -n $SGE_TASK_ID | \
    tail -n 1)
sample=$(basename $samplepath | awk '{ gsub("_WGS", "") ; print $0 }')

echo $sample

if [ -f $outdir/$sample.bam ]; then
   echo Sample $sample has been processed. Exiting.
   exit 0
fi

read1fqpaired=$scratchdir/$sample.paired.R1.fastq.gz
read2fqpaired=$scratchdir/$sample.paired.R2.fastq.gz

echo Aligning $sample with bowtie2
bowtie2 --seed 42 --very-fast --end-to-end --threads $CORES \
    -x $refgenome -1 $read1fqpaired -2 $read2fqpaired | \
    samtools view -1bS > $outdir/$sample.bam
