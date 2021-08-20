#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=100G
#$ -l mem_free=4.5G
#$ -l h_vmem=4.5G
#$ -l h_rt=96:00:00
#$ -pe local 8
#$ -t 1-41

## the following line is needed for sambamba
PATH="/users/scristia/.linuxbrew/bin:$PATH"
module load samtools
CORES=8

umask g+w

# Inputs
fqdir="../fastq"
bamdir="../bam"
beddir="../bed"
scratchdir="tmp"
## align with alternates
refgenome="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/genome/hg19/hg19"
outdirfastp=$bamdir/fastp

samplepath=$(find $fqdir -maxdepth 1 -name "*.fastq.gz" | \
    sed 's/.\{12\}$//' | \
    sort -u | \
    head -n $SGE_TASK_ID | \
    tail -n 1)
sample=$(basename $samplepath | awk '{ gsub("_WGS", "") ; print $0 }')

echo $sample

if [ -f $beddir/$sample.bed ]; then
   echo Sample $sample has been processed. Exiting.
   exit 0
fi

## Switch entirely to sambamba?
# sambamba view -h -S -f bam -l 1 - > $bamdir/$sample.bam

## Should we save the bam file with or without duplicates? -r option
## removes rather than marks duplicates. Would prefer to only keep one saved.
echo Flagging duplicates with SAMBAMBA.
sambamba markdup -t $CORES -l 1 --tmpdir=$TMPDIR \
    --hash-table-size=10000000 \
    --overflow-list-size=10000000 \
    --sort-buffer-size=2048 \
    -r $bamdir/$sample.bam $scratchdir/${sample}_rmdup.bam

## Sort with sambamba by read group then create bed file. 8th column is mapq
## score for read pair. Can filter with awk or downstream in R.
sambamba sort -n -t $CORES $scratchdir/${sample}_rmdup.bam -o /dev/stdout | \
    samtools fixmate -r - - |
    bamToBed -bedpe -i - | \
    awk '($2 != -1) && ($5 != -1) && ($1 == $4)' | \
    cut -f 1,2,6,8 | \
    sort -T $TMPDIR -k1,1 -k2,2n -S 5G > $beddir/$sample.bed

## Can pipe the above into bedtools nuc to get nucleotide frequencies for each
## base fragment. Takes a few hours and will be much faster if we use
## bedtools intersect -c compared to a sorted bed file of every G/C position.

#bedtools nuc -fi $hg19fa -bed - > $beddir/$sample.bed

## rm $scratchdir/${sample}_rmdup.bam $read1fqpaired $read2fqpaired

echo Done processing $sample!
