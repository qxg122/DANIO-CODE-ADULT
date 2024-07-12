#!/bin/bash
#SBATCH -J fastqdump
#SBATCH --array 0-4
#SBATCH --time=09:00:00

module purge; module load bluebear
module load SRA-Toolkit/2.10.9-gompi-2020b

while IFS= read -r line; do
  arr+=("$line")
done < ./Qianhong/SRR_Acc_List.txt

acc=${arr[$SLURM_ARRAY_TASK_ID]}

cd ./Qianhong/fastq

fasterq-dump ${acc} -t ./Qianhong/fastq
