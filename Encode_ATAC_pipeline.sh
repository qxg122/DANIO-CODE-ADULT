#!/bin/bash
#SBATCH -J atacpipeline
#SBATCH --array 0-12

cd ./Qianhong/

module purge; module load bluebear
module load Java/11.0.2

while IFS= read -r line; do
  arr+=("$line")
done < ./Qianhong/atac_sampleList.txt

atac=${arr[$SLURM_ARRAY_TASK_ID]}

caper hpc submit ./atac-seq-pipeline/atac.wdl -i ./Qianhong/atac_input/${atac} --singularity --leader-job-name r
