#!/usr/bin/env python3

import os, sys
import numpy as np
import scipy.stats

looplen = []
cosmut = []
totmut = []
unimut = []
clinmut = []
for line in open('../data/converge.tsv', 'r'):
    if 'GENE' not in line:
        start = int(line.split('\t')[22])
        end = int(line.split('\t')[23])
        tm = line.split('\t')[8]
        signal = int(line.split('\t')[9])
        signalp_TM = float(line.split('\t')[17])
        signalp_noTM = float(line.split('\t')[18])
        total_mutations = int(line.split('\t')[24])
        uniprot_mutations = int(line.split('\t')[25])
        clivar_mutations = int(line.split('\t')[26])
        cosmic_mutations = int(line.split('\t')[27])
        location = line.split('\t')[15]

        if location == 'Cytoplasmic' and tm == 'YES' and signal == 1 and signalp_TM <= 0.5 and signalp_noTM >= 0.5:
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
