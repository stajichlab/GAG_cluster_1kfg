#!/usr/bin/bash
#SBATCH -p intel --mem 64gb --out logs/clinker.log -N 1 -n 8 

module unload miniconda2
module load miniconda3
conda activate clinker

clinker -p cluster_contigs/CL1/* -p 1kfg_CL1_GAG.html
clinker -p cluster_contigs/CL2/* -p 1kfg_CL2_GAG.html

