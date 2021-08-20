#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=40G
#$ -l h_vmem=40G
#$ -l h_fsize=100G
#$ -l h_rt=24:00:00
#$ -M aannapr1@jhmi.edu
#$ -t 1-41



module load conda_R/4.0.x

fragdir="../granges"
outdir="../gc-counts"
mkdir -p $outdir

samplepath=$(find $fragdir -maxdepth 1 -name "*.rds"  | sort -u | head -n $SGE_TASK_ID | tail -n 1)
sample=$(basename ${samplepath//.rds})
out_file="${outdir}/${sample}_gc.csv"
if [ ! -f $out_file ]; then
    Rscript 02-gc_counts.r --id $samplepath --outdir $outdir
fi ;
