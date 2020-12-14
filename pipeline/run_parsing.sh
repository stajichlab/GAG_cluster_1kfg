#!/usr/bin/bash
#SBATCH -p short -N 1 -n 48 --mem 96gb --out logs/parsing.log
CPU=48
module unload perl
module load parallel
module load fasta

parallel -j $CPU ./scripts/build_cluster_genbank.py  --query {} ::: $(tail -n +2 jgi_names.tab | cut -f1)

parallel -j $CPU ./scripts/update_gbk_genes.py --query ::: $(ls cluster_contigs/CL[12]/*.gbk) 

mkdir -p cluster_by_subphyla/CL1 cluster_by_subphyla/CL2
./scripts/organize_by_clade.py cluster_contigs/CL1 cluster_by_subphyla/CL1
./scripts/organize_by_clade.py cluster_contigs/CL2 cluster_by_subphyla/CL2
