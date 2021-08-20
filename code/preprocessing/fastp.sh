#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=100G
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_rt=96:00:00
## #$ -pe local 8
#$ -t 1-41

module load bowtie2
module load samtools

PATH="/users/scristia/.linuxbrew/bin:$PATH"
FASTP=/users/scristia/fastp
CORES=1 ##CORES=8

umask g+w

# Inputs
fqdir="../fastq"
outdir="../bam"
beddir="../bed"
scratchdir="tmp"
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

read1fq=$samplepath\_R1.fastq.gz
read2fq=$samplepath\_R2.fastq.gz

# Trimming with fastp
echo Trimming $sample with fastp
read1fqpaired=$scratchdir/$sample.paired.R1.fastq.gz
read2fqpaired=$scratchdir/$sample.paired.R2.fastq.gz

if [ -f $read1fqpaired ]; then
   echo Sample $sample has been processed. Exiting.
   exit 0
fi

json=$outdirfastp/$sample.json
html=$outdirfastp/$sample.html
$FASTP -i $read1fq -I $read2fq -o $read1fqpaired -O $read2fqpaired \
    -Q -z 1 --detect_adapter_for_pe -j $json -h $html -w $CORES
