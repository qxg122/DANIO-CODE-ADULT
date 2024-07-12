#!/bin/bash
#SBATCH --mem=10G
#SBATCH -c 8
#SBATCH -o slurm-%J.out
#SBATCH --time=09:00:00

module purge; module load bluebear
module load Java/11.0.2

state_num=10
echo $state_num

while IFS= read -r line; do
  arr+=("$line")
done < ./Qianhong/tissue.txt

acc=${arr[$SLURM_ARRAY_TASK_ID]}

set -x

cd ./Qianhong/ChromHMM_file/${acc}

java -mx10G -jar ./Qianhong/ChromHMM/ChromHMM.jar BinarizeBam ./Qianhong/ChromHMM/CHROMSIZES/danRer11.txt . cell_mark_filetable.txt binary_bam_out

java -mx10G -jar ./Qianhong/ChromHMM/ChromHMM.jar LearnModel -p 10 -printposterior -printstatebyline binary_bam_out ChromHMM_out $state_num danRer11

set +x
