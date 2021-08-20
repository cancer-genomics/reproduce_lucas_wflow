#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=50G  ## this is high RAM. Observed 30.2G in AHN
#$ -l h_vmem=50G
#$ -l h_fsize=100G
#$ -l h_rt=24:00:00
#$ -M aannapr1@jhmi.edu
#$ -t 40-41


module load conda_R/4.0.x

CWD=$PWD
beddir="../bed"
fragdir="../granges"
mkdir -p $fragdir $beddir

cd $beddir
input=$(ls -1v *.bed | head -n $SGE_TASK_ID | tail -n 1)
sample=${input//.bed}

cd $CWD
bed_file=${beddir}/${input}
frag_file=${fragdir}/${sample}.rds
if [ ! -f $frag_file ]; then
 Rscript 01-bed_to_granges.r --id $bed_file --outdir $fragdir
 echo Finished creating fragments!
fi ;
