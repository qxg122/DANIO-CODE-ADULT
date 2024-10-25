#!/bin/bash
#SBATCH -J homer
#SBATCH --array 0-12
#SBATCH --mem=10G
#SBATCH -c 8
#SBATCH -o slurm-%J.out
#SBATCH --time=09:00:00

cd /rds/projects/m/muellerf-cage/Qianhong/Adult_Zebrafish_CREs_Annotation/

module purge; module load bluebear
module load Java/11.0.2

while IFS= read -r line; do
  arr+=("$line")
done < /rds/projects/m/muellerf-cage/Qianhong/Adult_Zebrafish_CREs_Annotation/tissues.txt

tissue=${arr[$SLURM_ARRAY_TASK_ID]}

findMotifsGenome.pl /rds/projects/m/muellerf-cage/Qianhong/Adult_Zebrafish_CREs_Annotation/${tissue}_PADREs_formal_idr_github_10_chrstart.bed.sorted.bed danRer11 ${tissue}_CRE_motif/ -size 200