#!/usr/bin/env python3

import os, re, csv, sys, argparse
from operator import itemgetter, attrgetter
from Bio import SeqIO
from Bio import bgzf

argparser = argparse.ArgumentParser(description='Parse one or many sequence files for GAG clusters.')
extmap = {'phmmer.domtbl': 'PHMMER',
        'TBLASTN.tab': 'TBLASTN',
        'TFASTX.tab':  'TFASTX' }
all = ['uge3','sph3','ega3','agd3','gtb3']
cluster_expected = { 'CL1': ['gtb3','agd3','ega3','sph3','uge3'], #,'uge3'
                     'CL2': ['gtb3','sph3','uge3'] } # ,'uge3'
argparser.add_argument('--indir',default='search',
                    help='input directory')
argparser.add_argument('--maxdist',default=15000,
                        help='maximum distance between genes to still be classified a cluster')
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
        rows = []
        for row in incsv:
            row[8] = int(row[8])
            row[9] = int(row[9])
            rows.append(row)

    last_end = None
    last_ctg = None
    cluster = []
    for row in sorted(rows,key=itemgetter(1,8)):
        q = row[0]
        h = row[1]
        if q not in phmmerhits:
            continue
        hstart = min(int(row[8]),int(row[9]))
        hend   = max(int(row[8]),int(row[9]))

        if not last_ctg or last_ctg != h:
            print("new cluster for {} {} {} {}".format(h,hstart,hend,q))
            cluster.append([h,hstart,hend,[q]])
        elif last_end and abs(hstart - last_end) > args.maxdist:
            print("new cluster [too far] for {} {} {} {}".format(h,hstart,hend,q))
            cluster.append([h,hstart,hend,[q]])
        else:
            cluster[-1][1] = min(cluster[-1][1],hstart)
            cluster[-1][2] = max(cluster[-1][2],hend)
            cluster[-1][3].append(q)

        last_end   = hend + 1
        last_ctg   = h

    n = 0
    for cl in cluster:
        cluster_classify = None
        ctgname = cl[0]

        for clname,genenames in cluster_expected.items():
            expected_len = len(genenames)
            gene_count   = []
            print("seen names is {}".format(cl[3]))
            for g in genenames:
                if g in cl[3]:
                    gene_count.append(g)
            if len(set(gene_count)) == expected_len:
                cluster_classify = clname
                print("gene_count is {} expected len is {}".format(gene_count,expected_len))
                break

        if not cluster_classify:
            print("--> no cluster on {}:{}".format(species,ctgname))
            continue

        outdir = os.path.join(args.clusterout,cluster_classify)
        if ctgname not in record_dict:
            print("cannot find {} in record for {}".format(ctgname,species))
        seq = record_dict[ctgname]
        ctgname = cl[0]
        left    = cl[1]
        right   = cl[2]
        # space out a bit more
        left -= int(args.maxdist / 2)
        if left < 1:
            left = 1
        right += int(args.maxdist / 2)
        if right > len(seq):
            right = len(seq)

        print("slicing {}_{}:{}..{}".format(species,ctgname,left,right))
        slice = seq[ left:right ] # cut a slice out
        slice.id = "{}_{}".format(species,seq.id)
#        for feature in slice.features:
#            if "gene" in feature.qualifiers:
#                lt =  features.qualifiers.get("locus_tag")
#                gene = feature.qualifiers.get("gene")
#                feature.qualifiers["locus_tag"] = gene
#                feature.qualifiers["gene"] = lt
        SeqIO.write(slice, os.path.join(args.clusterout,
            cluster_classify,"{}.{}_cls{}.gbk".format(species,ctgname,n)), "genbank")
        n += 1
