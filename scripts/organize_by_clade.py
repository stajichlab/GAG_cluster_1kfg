#!/usr/bin/env python3
import csv, sys, re, os

if len(sys.argv) < 3:
    print("Expected 2 arguments {}".format(sys.argv))
    sys.exit(0)

indir = sys.argv[1]
outdir = sys.argv[2]

if not os.path.isdir(outdir):
    os.mkdir(outdir)

taxonomytable="jgi_names_taxonomy.csv"
taxlevel_split = "class"
clades = {}
taxon2clade = {}
with open(taxonomytable,"rt") as fh:
    csvin = csv.reader(fh)
    hdr = next(csvin)
    for row in csvin:
        prefix=row[0]
        species=row[1]
        straininfo=row[2]
        taxonomy = row[3]
        clade = ""
        if ( taxonomy.startswith("phylum:Microsporidia") or taxonomy.startswith("phylum:Cryptomycota")):
            # skip Microsporidia
            continue
        for taxnode in taxonomy.split(";"):
            if taxnode.startswith(taxlevel_split):
                tn = taxnode.split(":")
                clade = tn[1]
                break
            elif taxnode.startswith("family:"):
                tn = taxnode.split(":")
                clade = tn[1]
                break

        taxon2clade[prefix] = [clade, species, straininfo]

        if clade not in clades:
            clades[clade] = []

        clades[clade].append([prefix,species,straininfo])

nameparse = re.compile(r'^([^\.]+)\.(\S+)')
MASTERSRC = ""
MASTERNAME=""
for fname in os.listdir(indir):
    if not fname.endswith(".gbk"):
        continue

    m = nameparse.match(fname)
    prefix = ""
    rest   = ""
    if m:
        prefix = m.group(1)
        rest   = m.group(2)
    else:
        print("cannot see a prefix in {}".format(fname))
        continue

    if prefix in taxon2clade:
        clade = taxon2clade[prefix][0]
        odir = os.path.join(outdir,clade)
        if not os.path.isdir(odir):
            os.mkdir(odir)
        species = taxon2clade[prefix][1]
        strain  = taxon2clade[prefix][2]
        strain  = re.sub(r'\/','_',strain)
        fullname = re.sub(r'\s+','_',species + " " + strain) + "." + rest
        sourcepath = os.path.realpath(os.path.join(indir,fname))
        destpath   = os.path.join(odir,fullname)
        print("ln -s {} {}".format(sourcepath, destpath))
        if not os.path.islink(destpath):
            os.symlink( sourcepath, destpath )
        if species == "Aspergillus fumigatus":
            MASTERSRC=sourcepath
            MASTERNAME=fullname
for clade in clades:
    if not os.path.isdir(os.path.join(outdir,clade)):
        continue
    dest = os.path.join(outdir,clade,MASTERNAME)
    if not os.path.islink(dest):
        print("Adding MASTER {} -> {}".format(MASTERSRC,dest))
        os.symlink(MASTERSRC,dest)
