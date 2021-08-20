#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=1G
#$ -l h_vmem=1G
# Job resource option: max runtime
#$ -l h_rt=96:00:00
qsub fastp.sh
qsub -hold_jid_ad fastp.sh align.sh
qsub -hold_jid_ad align.sh post_alignment.sh
qsub -hold_jid_ad post_alignment.sh,align.sh bed_to_granges.sh
qsub -hold_jid_ad post_alignment.sh,bed_to_granges.sh,align.sh gc_counts.sh
qsub -hold_jid post_alignment.sh,bed_to_granges.sh,gc_counts.sh,align.sh bin_corrected.sh
