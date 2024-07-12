#!/bin/bash
#SBATCH -J chipipeline
#SBATCH --array 0-3

cd ./Qianhong/

module purge; module load bluebear
module load Java/11.0.2

while IFS= read -r line; do
  arr+=("$line")
done < ./Qianhong/chip_sampleList.txt

chip=${arr[$SLURM_ARRAY_TASK_ID]}

caper hpc submit ./chip-seq-pipeline2/chip.wdl -i ./Qianhong/chip_input/${chip} --singularity --leader-job-name v
