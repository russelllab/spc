#!/usr/bin/env python3

## Script to plot SPC
import os, sys, gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pylab
import numpy as np

class spc:
    def __init__(self, id):
        self.id = ''
        self.acc = '-'
        self.gene = '-'
        self.kw = 'noTM'
        self.signal = 0
        self.topo = {}
        self.fasta = ''
        self.type = 'NA'
        self.positions = {}
        self.transmem = {}
        self.pos_tm = 'NA'
        self.pos_notm = 'NA'
        self.prob_tm = 0
        self.prob_notm = 0
        self.pred_tm = ''
        self.pred_notm = ''
        self.mut = []
        self.uniprot_mut = []
        self.cosmic_mut = []
        self.clinvar_mut = []
        self.topfind = []

x = []
flag = 0
dic = {}
for line in gzip.open('../../DB/uniprot/uniprot_sprot_human.dat.gz', 'rt'):
    #print (line)
    if line[:2] == 'ID':
        id = line.split()[1]
        dic[id] = spc(id)
        type = 'Unknown'
        flag = 0
        count = 0
        sq = 0
    elif line[:2] == 'AC' and dic[id].acc == '-':
        dic[id].acc = line.split('AC')[1].replace(' ', '').replace('\n', '').split(';')[0]
    elif line[:2] == 'GN' and dic[id].gene == '-':
        #print (line)
        dic[id].gene = line.split('=')[1].split(' ')[0].replace(' ', '').replace(',', '').replace(';', '').replace('\n', '')
    elif line[:2] == 'SQ':
        sq = 1
        continue
    elif line[:2] == '//':
        sq = 0
        continue
    elif line[:2] == 'KW' and dic[id].kw == 'noTM':
        if 'Transmembrane' in line:
            dic[id].kw = 'TM'
    elif line[:2] == 'FT':
        if line.split()[1] in ['SIGNAL']:
            dic[id].signal = 1
        if line.split()[1] in ['SIGNAL', 'TRANSMEM', 'TOPO_DOM']:
            topo = line.split()[1]
            if topo == 'TOPO_DOM':
                count += 1
                start = line.split()[2].split('.')[0]
                end = line.split()[2].split('.')[-1].replace('\n', '')
                if '?' not in start and '?' not in end:
                    start = int(start)
                    end = int(end)
                    flag = 1
                    #if id == 'CXB1_HUMAN':
                    #print (count, id, count, topo)
                    dic[id].topo[topo+str(count)] = [start, end]
                    #dic[id].topo[topo+str(count)].append(start)
                    #dic[id].topo[topo+str(count)].append(end)
        else:
            if '/note=' in line.split()[1] and flag == 1:
                type = line.split()[1].split('"')[1]
                dic[id].topo[topo+str(count)].append(type)
                flag = 0

        #break
    if sq == 1:
        dic[id].fasta += line.replace('\n', '').replace(' ', '')

print (dic['CXB1_HUMAN'].topo)

map = {}
mapgid = {}
for line in gzip.open('/home/gurdeep/projects/DB/uniprot/HUMAN_9606_idmapping.dat.gz', 'rt'):
    #print (line)
    if line.split('\t')[1] == 'UniProtKB-ID':
        acc = line.split('\t')[0]
        id = line.split('\t')[2].replace('\n', '')
        map[acc] = id
    elif line.split('\t')[1] == 'Gene_Name':
        if line.split('\t')[0] == acc:
            gene = line.split('\t')[2].replace('\n', '')
            if gene not in mapgid:
                mapgid[gene] = id

print (mapgid['RHBDF2'])

'''
for line in open('../data/top_find_output.tsv', 'rt'):
    acc = line.split(' ')[0].replace(' ', '')
    position = line.split(',')[2].replace(' ', '')
    id = map[acc]
    if position not in dic[id].topfind:
        dic[id].topfind.append(position)
'''

## COSMIC mutations
for line in gzip.open('../data/cosmic.tsv.gz', 'rt'):
    if line[0] != '#':
        gene = line.split('\t')[0]
        if gene in mapgid:
            id = mapgid[gene]
            if id in dic:
                #dic[id].cosmic_mut = int(line.split('\t')[1].replace('\n', ''))
                num_sampleid = int(line.split('\t')[2].replace('\n', ''))
                if num_sampleid >= 3:
                    dic[id].cosmic_mut.append(line.split('\t')[1].split('.')[1])
                    dic[id].mut.append(line.split('\t')[1].split('.')[1])

## ClinVar mutations
for line in gzip.open('../data/clinvar_missense.tsv.gz', 'rt'):
    if line[0] != '#':
        gene = line.split('\t')[0]
        if gene in mapgid:
            id = mapgid[gene]
            if id in dic:
                dic[id].clinvar_mut.append(line.split('\t')[1].replace('\n', ''))
                dic[id].mut.append(line.split('\t')[1].replace('\n', ''))

## UniProt mutations
for line in gzip.open('../data/uniprot_2020-05_features_Hsa.tsv.gz', 'rt'):
    if line.split('\t')[1] == 'VARIANT' and '->' in line.split('\t')[5]:
        acc = line.split('\t')[0]
        #id = line.split('\t')[1]
        id = map[acc]
        if id in dic:
            #mut = line.split('\t')[0].split('/')[1].split()[0]
            #print (line)
            mutation = line.split('\t')[5].split('-')[0]
            mutation += str(line.split('\t')[4])
            mutation += line.split('\t')[5].split('>')[1]
            if 'in'  in line:
                dic[id].uniprot_mut.append(mutation)
                dic[id].mut.append(mutation)


for id in dic:
    dic[id].mut = list(set(dic[id].mut))

print (dic['CXB1_HUMAN'].mut)
#sys.exit()

def count_mutations(id, start, end, x):
    count = 0
    for mutation in x:
        #print (id, mutation, mutation[1:-1])
        position = ''
        for char in mutation:
            if char.isnumeric() == True:
                position += char
        position = int(position)
        if position >= start and position <= end:
            count += 1

    return (count)

df = pd.DataFrame()
data = []
#l = 'ID\tGENE\tACC\tTM protein\tNum MUT\tPosition\tSignalP_TM\tSignalP_noTM\tN-terminus\tTopology\n'
l = 'ID\tGENE\tACC\tTM protein\tUniProt MUT\tClinVar MUT\tCOSMIC MUT\tPosition\tSignalP_TM\tSignalP_noTM\tN-terminus\tTopology\n'
total_count = 0
count_excl = 0
count_type1 = 0
count_type2 = 0
count_type3 = 0
count_dom1 = 0
for id in dic:
    if dic[id].signal in [0, 1]:
    #if dic[id].signal in [0] and dic[id].kw == 'YES':
    #if dic[id].signal in [0] and dic[id].type == 'Cytoplasmic':
    #if dic[id].signal in [0]:
        #print (id, dic[id].acc)
        for topo in dic[id].topo:
            start = dic[id].topo[topo][0]
            end = dic[id].topo[topo][1]
            type = dic[id].topo[topo][2]
            name = id+'_'+str(topo)

            #print (name)
            '''
            if start <= 5:
                start_seq = 0
                end_seq = 70
            else:
                start_seq = start-5-1
                end_seq = start_seq+70
                if end_seq > len(dic[id].fasta):
                    start_seq = len(dic[id].fasta) - 70
                    end_seq = len(dic[id].fasta)

                    check = []
                    for tm_check in dic[id].transmem:
                        type_check, start_check, end_check = dic[id].transmem[tm_check]
                        #print (type_check, start_check, end_check, start_seq, end_seq, len(dic[id].fasta))
                        if (start_check<=end_seq and start_check>=start_seq) or (end_check<=end_seq and end_check>=start_seq):
                            row =[]
                            #print (type_check, start_check, end_check, end='\t')
                            row.append(type_check)
                            row.append(start_check)
                            row.append(end_check)
                            check.append(row)
                    #if len(check) == 1 and type != 'Unknown':
                    if len(check) > 1:
                        count_excl += 1
                        print (id, check)
                        continue
                        sys.exit()

            '''
            '''
            prob_tm = 0.0
            prob_notm = 0.0
            if os.path.isfile('../data/part_B_signalp/'+name+'_long_noTM.out'):
                #print ('YES')
                for line in open('../data/part_B_signalp/'+name+'_long_noTM.out', 'r'):
                    if line.split()[0] == 'max.' and line.split()[1]=='Y':
                        prob_notm = float(line.replace('\n','').split()[3])
                        pos_notm = int(line.replace('\n','').split()[2]) + map_start[name] - 1
                        dic[id].pos_notm = pos_notm
                    if line.split()[0] == 'D':
                        pred_notm = line.replace('\n','').split()[-1]
                for line in open('../data/part_B_signalp/'+name+'_long_TM.out', 'r'):
                    if line.split()[0] == 'max.' and line.split()[1]=='Y':
                        prob_tm = float(line.replace('\n','').split()[3])
                        pos_tm = int(line.replace('\n','').split()[2]) + map_start[name] - 1
                        dic[id].pos_tm = pos_tm
                    if line.split()[0] == 'D':
                        pred_tm = line.replace('\n','').split()[-1]
            else:
                print ('No', '../data/part_B_signalp/'+name+'_long_noTM.out')
                #sys.exit()
            '''

            row = []
            #print (name, dic[id].gene, dic[id].acc, dic[id].kw, dic[id].mut, prob_tm, prob_notm, dic[id].type, dic[id].topo)
            #l += id+'\t'+ dic[id].gene+'\t'+ dic[id].acc+'\t'+ dic[id].kw+'\t'+ str(dic[id].mut)+'\t'+ str(dic[id].pos_notm) + '\t' + str(dic[id].prob_tm)+'\t'+ str(dic[id].prob_notm)+'\t'+ dic[id].type+'\t'+ dic[id].topo+'\n'
            row.append(id)
            row.append(dic[id].gene)
            row.append(dic[id].acc)
            row.append(name)
            row.append(topo)
            row.append(start)
            row.append(end)
            #row.append(dic[id].pos_notm)
            row.append(dic[id].kw)
            row.append(dic[id].signal)
            row.append(count_mutations(id, int(start), int(end), dic[id].mut))
            row.append(count_mutations(id, int(start), int(end), dic[id].uniprot_mut))
            row.append(count_mutations(id, int(start), int(end), dic[id].clinvar_mut))
            row.append(count_mutations(id, int(start), int(end), dic[id].cosmic_mut))
            #row.append(5)
            #row.append(type)

            if type in ['Cytoplasmic', 'Extracellular', 'Lumenal', 'Nuclear']:
                row.append(type)
            else:
                row.append('Other')

            data.append(row)


df = pd.DataFrame(data, columns=['ID', 'GENE', 'ACC', 'Name', 'TOPO', 'Start', 'End', 'Transmembrane','SIGNALp', 'Total_MUT', 'Uniprot_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'N-terminus'])
df.to_csv('../data/loop_length_mutations.tsv', sep='\t')
print ('completed')
