#!/usr/bin/env python3

import os, re, csv, sys, argparse
from Bio import SeqIO
from Bio import bgzf

argparser = argparse.ArgumentParser(description='Parse one or many sequence files for GAG clusters.')
extmap = {'phmmer.domtbl': 'PHMMER',
        'TBLASTN.tab': 'TBLASTN',
        'TFASTX.tab':  'TFASTX' }
all = ['uge3','sph3','ega3','agd3','gtb3']
cluster_expected = { 'CL1': ['gtb3','agd3','ega3','sph3','uge3'],
                     'CL2': ['gtb3','sph3','uge3'] }
argparser.add_argument('--indir',default='search',
                    help='input directory')
argparser.add_argument('--lookupname',default='jgi_names.tab',
                        help='Species names lookup TSV file')
argparser.add_argument('--gbkdir',default='GBK',
                        help='GBK (bgzip files) directory')
argparser.add_argument('--clusterout',default='cluster_contigs',
                        help='Out directyory for writing cluster contigs in gbk')
argparser.add_argument('--query','--species',help='Name of one species file to process')

args = argparser.parse_args()

if not os.path.isdir(args.clusterout):
    os.mkdir(args.clusterout)
for clname in cluster_expected:
    cldir = os.path.join(args.clusterout,clname)
    if not os.path.isdir(cldir):
        os.mkdir( cldir )

hits = {}
pref2species = {}
sppat = re.compile(r'^(\S+)\.GAGcluster')

with open(args.lookupname,"rt") as fh:
    rdr = csv.reader(fh,delimiter="\t")
    header = next(rdr)
    for row in rdr:
        pref2species[row[0]] = row[1]


for fname in os.listdir(args.indir):
    m = sppat.search(fname)
    species = ""
    if m:
        species = m.group(1)
    else:
        print("cannot match {}".format(fname))
        continue
    if args.query and species != args.query:
#        print("searching for one species {}".format(args.query))
        continue

    nameext = ".".join(fname.split(".")[-2:])
    if nameext in extmap:
        type = extmap[nameext]
        #print(extmap[nameext],nameext,fname)
        if species in hits:
            hits[species][type] = os.path.join(args.indir,fname)
        else:
            hits[species] = {type: os.path.join(args.indir,fname)}

for species in hits:
    cluster = {}
    record_dict = SeqIO.index(os.path.join(args.gbkdir,"{}.gbf.gz".format(species)),'genbank')
    phmmerhits = {}
    # is there an agd3 or gtb3 in the protein search?
    with open(hits[species]['PHMMER']) as infh:
            for line in infh:
                row = line.split()
                if row[0].startswith("#"):
                    continue
                hit = row[0]
                query = row[3]
                phmmerhits[query] = hit

    with open(hits[species]['TFASTX']) as fh:
        incsv = csv.reader(fh, delimiter="\t")
        seen = {}
        ctgrange = {}
        for row in incsv:
            q = row[0]
            h = row[1]
            if q not in phmmerhits:
                continue
            hstart = min(int(row[8]),int(row[9]))
            hend   = max(int(row[8]),int(row[9]))

            if h not in cluster:
                cluster[h] = []

            cluster[h].append([q,int(row[8]),int(row[9])])

            if h not in seen:
                seen[h]     = {q: 1}
                ctgrange[h] = [hstart,hend]
            else:
                seen[h][q]     = 1
                ctgrange[h][0] = min(hstart,ctgrange[h][0])
                ctgrange[h][1] = max(hend,ctgrange[h][1])

        for contig in seen:
            cluster_classify = None
            for clname,genenames in cluster_expected.items():
                expected_len = len(genenames)
                gene_count   = 0
                for g in genenames:
                    if g in seen[contig]:
                        gene_count += 1
                if  gene_count == expected_len:
                    cluster_classify = clname
                    break

            if not cluster_classify:
                print("no cluster collection of genes found in {}".format(contig))
                continue

            outdir = os.path.join(args.clusterout,cluster_classify)
            seq = record_dict[contig]

            print("seq is {}".format(seq.id))
            print("slicing for {}_{}:{}..{}".format(species,contig,ctgrange[contig][0],ctgrange[contig][1]))
            slice = seq[ ctgrange[ contig][0]:ctgrange[contig][1] ] # cut a slice out
            slice.id = "{}_{}".format(species,seq.id)
            SeqIO.write(slice, os.path.join(args.clusterout,cluster_classify,"{}_{}.gbk".format(species,contig)), "genbank")

#            if 'agd3' in seen[contig] and 'gtb3' in seen[contig]:
#
#            if 'gtb3' in seen[contig] and len(seen[contig].keys()) > 1:
#                seq = record_dict[contig]
#                print("seq is {}".format(seq.id))
#                print("slicing for {}_{}:{}..{}".format(species,contig,ctgrange[contig][0],ctgrange[contig][1]))
#                slice = seq[ ctgrange[contig][0]:ctgrange[contig][1] ] # cut a slice out
#                slice.id = "{}_{}".format(species,seq.id)
#                SeqIO.write(slice, os.path.join(args.clusterout,"{}_{}.gbk".format(species,contig)), "genbank")
