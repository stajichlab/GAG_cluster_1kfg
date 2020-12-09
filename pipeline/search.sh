#!/usr/bin/bash
#SBATCH -p short -N 1 -n 4 --mem 2gb --out logs/search.%a.log

module load ncbi-blast/2.9.0+
module load fasta
module load hmmer/3
SAMPFILE=jgi_names.tab
QUERY=query/GAGcluster.fas
QUERYHARD=query/GAGcluster.hard.fas
QBASE=$(basename $QUERY .fas)
DB=DNA
PEPDB=PEP
EVALUE=1e-4
EVALUE_PHMMER=1e-5
OUT=search
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPFILE"
  exit
fi

tail -n +2 $SAMPFILE | sed -n ${N}p | while read BASE SPECIES
do
	if [ ! -s $OUT/${BASE}.$QBASE.TBLASTN.tab ]; then
		tblastn -query $QUERY -db $DB/${BASE}.nt.fasta -num_threads $CPU -evalue $EVALUE -outfmt 6 -out $OUT/${BASE}.$QBASE.TBLASTN.tab
	fi
	if [ ! -f $OUT/${BASE}.$QBASE.TFASTX.tab ]; then
		tfasty36 -k 1000 -S -T $CPU -E $EVALUE -m 8c $QUERYHARD $DB/${BASE}.nt.fasta > $OUT/${BASE}.$QBASE.TFASTX.tab 
	fi
	if [ ! -f $OUT/${BASE}.$QBASE.phmmer.domtbl ]; then
		phmmer --cpu $CPU -E $EVALUE_PHMMER --domtbl $OUT/${BASE}.$QBASE.phmmer.domtbl $QUERY $PEPDB/${BASE}.aa.fasta | gzip -c > $OUT/${BASE}.$QBASE.phmmer.out.gz
	fi
done
