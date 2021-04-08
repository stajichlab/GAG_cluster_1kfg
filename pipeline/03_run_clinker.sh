#!/usr/bin/bash
#SBATCH -p intel --mem 64gb --out logs/clinker.log -N 1 -n 5

module unload miniconda2
module load miniconda3
conda activate clinker
module load parallel
which perl
PREF=1kfg
TOP=cluster_by_subphyla
# for CL in CL1 CL2
# do
#   echo "$CL"
#   for d in $(ls $TOP/$CL | grep -v \.gbk)
#   do
#     if [ -d $TOP/$CL/$d ]; then
#       echo "$CL/$d"
#       clinker -p $TOP/$CL/$d/*.gbk -p ${PREF}_${CL}_${d}.html
#     fi
#   done
# done

for CL in CL1_CL2 CL1 CL2
do
  parallel -j 5 clinker -p $TOP/$CL/{}/*.gbk -p ${PREF}_${CL}_{}.html ::: $(ls $TOP/$CL | grep -v \.gbk)
done
