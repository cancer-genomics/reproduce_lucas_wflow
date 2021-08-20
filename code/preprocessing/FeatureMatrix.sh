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
Rscript 06-FeatureMatrix.r LUCAS_Hiseq_55/zscores_hiseq55.csv LUCAS_Hiseq_55/coverage_hiseq55.csv LUCAS_Hiseq_55/allfeatures.csv
