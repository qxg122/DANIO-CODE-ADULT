#!/bin/bash
#SBATCH --mem=10G
#SBATCH -c 8
#SBATCH -o slurm-%J.out
#SBATCH --time=09:00:00
# Please state the number of array iteration at sbatch submission

module purge; module load bluebear
module load Java/11.0.2

# careful, this is hardcoded, depends on the order of stage input!!!!

#num_states=( 6 3 5 3 3 6 6 )
state_num=10
echo $state_num

cd /rds/projects/m/muellerf-cage/Qianhong/ChromHMM_file/chromHMM_combined/

set -x
java -mx10G -jar /rds/projects/m/muellerf-cage/Qianhong/ChromHMM/ChromHMM.jar BinarizeBam /rds/projects/m/muellerf-cage/Qianhong/ChromHMM/CHROMSIZES/danRer11.txt . combined_cell_mark_filetable.txt binary_bam_combined_out_10

java -mx10G -jar /rds/projects/m/muellerf-cage/Qianhong/ChromHMM/ChromHMM.jar LearnModel -p 10 -printposterior -printstatebyline binary_bam_combined_out_10 ChromHMM__combinedout_10 $state_num danRer11

set +x