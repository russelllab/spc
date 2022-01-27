#!/usr/bin/env python3

import os, sys
import numpy as np
import scipy.stats

looplen = []
cosmut = []
totmut = []
unimut = []
clinmut = []
for line in open('../data/loop_length_mutations.tsv', 'r'):
    if 'GENE' not in line:
        start = int(line.split('\t')[6])
        end = int(line.split('\t')[7])
        tm = line.split('\t')[8]
        signal = int(line.split('\t')[9])
        total_mutations = int(line.split('\t')[10])
        uniprot_mutations = int(line.split('\t')[11])
        clivar_mutations = int(line.split('\t')[12])
        cosmic_mutations = int(line.split('\t')[13])
        location = line.split('\t')[-1].replace('\n','')

        if location == 'Cytoplasmic' and tm == 'TM':
        #if location == 'Extracellular':
        #if True:
            looplen.append(end - start + 1)
            cosmut.append(cosmic_mutations)
            clinmut.append(clivar_mutations)
            unimut.append(uniprot_mutations)
            totmut.append(total_mutations)

print ('COSMIC', scipy.stats.pearsonr(looplen,cosmut))
print ('ClinVar', scipy.stats.pearsonr(looplen,clinmut))
print ('UniProt', scipy.stats.pearsonr(looplen,unimut))
print ('Total', scipy.stats.pearsonr(looplen,totmut))
