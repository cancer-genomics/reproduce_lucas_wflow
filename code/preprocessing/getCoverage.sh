#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=100G
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l h_rt=96:00:00
# Log output
#$ -o /dcs04/scharpf/data/annapragada/DELFI_pipeline_updates/logs

module load conda_R/4.0.x
Rscript 05-getCoverage.r LUCAS_Hiseq_55/coverage_hiseq55.csv /dcs04/scharpf/data/annapragada/DELFI_pipeline_updates/LUCAS_Hiseq_55/5mb_bins

