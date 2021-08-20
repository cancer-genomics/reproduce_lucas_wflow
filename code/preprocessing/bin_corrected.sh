#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -l h_fsize=100G
#$ -l h_rt=24:00:00
#$ -M aannapr1@jhmi.edu
#$ -t 345
# Log output
#$ -o /dcs04/scharpf/data/annapragada/DELFI_pipeline_updates/logs_55


rlib="${HOME}/Library/R/3.12-bioc-release-conda"
module load conda_R/4.0.x

CWD=$PWD
fragdir="/dcl02/leased/cglab/rscharpf/cristiano/projects/lucas/granges"
bindir="LUCAS_Hiseq_55/5mb_bins"
binfile="/dcl02/leased/cglab/rscharpf/cristiano/projects/lucas/bins/reference/bins_5mb.csv"
target="hi"
mkdir -p $bindir

samplepath=$(find $fragdir -maxdepth 1 -name "*.rds"  | sort -u | head -n $SGE_TASK_ID | tail -n 1)
sample=$(basename $samplepath | awk '{ gsub(".rds", "") ; print $0}')

frag_file=${fragdir}/${sample}.rds
bin_file=${bindir}/${sample}_5mb.csv
if [ ! -f $bin_file ]; then
    R_LIBS_USER=$rlib Rscript 03-bin_corrected.r --id $frag_file --outdir $bindir --bins $binfile --target $target
fi ;
