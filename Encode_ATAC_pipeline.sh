#!/bin/bash
#SBATCH -J atacpipeline
#SBATCH --array 0-12

cd /rds/projects/m/muellerf-cage/Qianhong/

module purge; module load bluebear
module load Java/11.0.2

while IFS= read -r line; do
  arr+=("$line")
done < /rds/projects/m/muellerf-cage/Qianhong/atac_sampleList.txt

atac=${arr[$SLURM_ARRAY_TASK_ID]}

caper hpc submit /rds/homes/q/qxg122/pipelines/atac-seq-pipeline/atac.wdl -i /rds/projects/m/muellerf-cage/Qianhong/atac_input/${atac} --singularity --leader-job-name r