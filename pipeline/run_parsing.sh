#!/usr/bin/bash
for n in $(tail -n +2 jgi_names.tab | cut -f1); do sbatch -p short -N 1 -n 1 --mem 2gb --wrap "./scripts/build_cluster_genbank.py  --query $n"; done
