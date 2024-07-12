#!/bin/bash
#SBATCH -J chipipeline
#SBATCH --array 0-3

cd /rds/projects/m/muellerf-cage/Qianhong/

module purge; module load bluebear
module load Java/11.0.2

while IFS= read -r line; do
  arr+=("$line")
done < /rds/projects/m/muellerf-cage/Qianhong/chip_sampleList_k9_me3.txt

chip=${arr[$SLURM_ARRAY_TASK_ID]}

caper hpc submit /rds/homes/q/qxg122/pipelines/chip-seq-pipeline2/chip.wdl -i /rds/projects/m/muellerf-cage/Qianhong/chip_input_k9/${chip} --singularity --leader-job-name v