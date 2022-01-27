#!/usr/bin/env python3

import os, sys, gzip
from Bio.PDB.Polypeptide import three_to_one

## ClinVar files
dic = {}
for line in gzip.open('../data/variant_summary.txt.gz', 'rt'):
    if line[0] != '#':
        gene = line.split('\t')[4]
        clinsig = line.split('\t')[6]
        mutation = line.split('\t')[2]
        type = line.split('\t')[1]
        if 'germline' in line and type == 'single nucleotide variant':
            if clinsig in ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']:
                if 'p.' in mutation:
                    mutation = mutation.split('p.')[1].split(')')[0]
                    if 'Ter' not in mutation and '=' not in mutation and 'X' not in mutation and '*' not in mutation:
                        if gene not in dic:
                            dic[gene] = []
                        #print (mutation[:3], mutation[3:-3], mutation[-3:], gene)
                        mutation =  three_to_one(mutation[:3].upper()) + mutation[3:-3] +  three_to_one(mutation[-3:].upper())
                        dic[gene].append(mutation)
                        dic[gene] = list(set(dic[gene]))

l = '#Gene\tMUT\n'
for gene in dic:
    for mutation in dic[gene]:
        l += gene + '\t' + mutation + '\n'

gzip.open('../data/clinvar_missense.tsv.gz', 'wt').write(l)
