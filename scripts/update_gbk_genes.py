#!/usr/bin/env python3
import os
import sys
import csv
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

searchcmd = "ssearch36 -k 1000 -m 8c -E 1e-3 -T 2 %s %s"

argparser = argparse.ArgumentParser(description='Rename genes in Genbank Contigs file.')
argparser.add_argument('--query', '--in',required=True,
                       help='Name of one species file to process')

argparser.add_argument('--db', '--querydb',required=False, default='query/GAGcluster.fas',
                       help='Query database for searches')

argparser.add_argument('--tmpdir',required=False,default="tmp_search",
                       help='Temporary Pep Folder')

args = argparser.parse_args()

if not os.path.isdir(args.tmpdir):
    os.mkdir(args.tmpdir)

if args.query.endswith(".gbk"):
    prefix = os.path.splitext(os.path.basename(args.query))[0]
    seqs = []
    for seq_record in SeqIO.parse(args.query, "genbank"):
        genes = {}
        pepsout = []
        for feature in seq_record.features:
            if feature.type in ["gene", 'mRNA', 'CDS']:
                locus = feature.qualifiers['locus_tag'][0]
                genename = feature.qualifiers['gene'][0]
#                        print("type is {} locus is {} gene is {}".format(feature.type,locus,genename))
                if feature.type == 'CDS':
                    pep = feature.qualifiers['translation'][0]
                    pepsout.append(
                        SeqRecord(Seq(pep), id=genename, description=locus))
                if genename not in genes:
                    genes[genename] = [locus]
#                        print("storing {} -> {}".format(genename,locus))
        pepfile = os.path.join(args.tmpdir, "%s.faa" % (prefix))
        SeqIO.write(pepsout, pepfile, "fasta")
        with os.popen(searchcmd % (pepfile, args.db)) as stream:
            tabparse = csv.reader(stream, delimiter="\t")
            for row in tabparse:
                q = row[0]
                hname = row[1]
                evalue = row[-2]
#                        print("q={} h={} e={}".format(q,hname,evalue))
#                        if q not in genes:
#                            print(" no record for {} in genes".format(q))
#                        else:
#                            print("storing {} -> {} {} - old value is {}".format(q,hname,evalue,genes[q]))
                genes[q].append(hname)
                genes[q].append(evalue)
#                        print("new value for {} -> {}".format(q,genes[q]))
        for feature in seq_record.features:
            if feature.type in ["gene", 'mRNA', 'CDS']:
                locusin = feature.qualifiers['locus_tag'][0]
                genename = feature.qualifiers['gene'][0]
                if genename in genes and len(genes[genename]) > 1:
                    feature.qualifiers['gene'] = genes[genename][1]

                feature.qualifiers['locus_tag'] = genename
        seqs.append(seq_record)
    os.rename(args.query, args.query + ".old")
    SeqIO.write(seqs, args.query, "genbank")
