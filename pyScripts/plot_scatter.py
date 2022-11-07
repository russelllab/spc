## Script to generate and plot N-terminal
## analysis data
import os, sys, gzip, re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pylab
import numpy as np


def fetchEnzymeCleavageSites(fasta, spcCleavageSite, enzyme):
    #pattern = motif[identifier].pattern
    print (fasta)
    #pattern = '[R]'
    #trypsin = '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
    pattern = enzyme
    x = re.finditer(pattern, fasta, flags=0)
    s = None; e = None; flag = 0
    for val in x:
        s, e = val.span()
        s += 1
        e += 1
        if int(spcCleavageSite) < int(e):
            flag = 1
            break
    if flag == 0:
        e = None
    print(e)
    return (e)
    #print (x)
    #sys.exit()

class spc:
    def __init__(self, id):
        self.id = ''
        self.acc = '-'
        self.gene = '-'
        self.kw = 'NO'
        self.signal = 0
        self.topo = ''
        self.fasta = ''
        self.type = '-'
        self.uniprotkb = '-'
        self.pos_tm = 0
        self.pos_notm = 0
        self.prob_tm = 0
        self.prob_notm = 0
        self.pred_tm = ''
        self.pred_notm = ''
        self.mut = []
        self.uniprot_mut = []
        self.cosmic_mut = []
        self.clinvar_mut = []
        #self.clinvar_mut = 0
        self.topfind = []

## Select the species
#species = 'caeel'
species = 'human'
sprot = 1
trembl = 0

if sprot == 1 and trembl == 1:
    if species == 'human':
        files = ['../../DB/uniprot/uniprot_sprot_human.dat.gz', '../../DB/uniprot/uniprot_trembl_human.dat.gz']
    else:
        files = ['../../DB/uniprot/uniprot_sprot_caeel.dat.gz', '../../DB/uniprot/uniprot_trembl_caeel.dat.gz']
    fastaFiles = ['../../DB/uniprot/uniprot_sprot.fasta.gz', '../../DB/uniprot/uniprot_trembl.fasta.gz']
elif sprot == 1 and trembl == 0:
    if species == 'human':
        files = ['../../DB/uniprot/uniprot_sprot_human.dat.gz']
    else:
        files = ['../../DB/uniprot/uniprot_sprot_caeel.dat.gz']
    fastaFiles = ['../../DB/uniprot/uniprot_sprot.fasta.gz']
elif sprot == 0 and trembl == 1:
    if species == 'human':
        files = ['../../DB/uniprot/uniprot_trembl_human.dat.gz']
    else:
        files = ['../../DB/uniprot/uniprot_trembl_caeel.dat.gz']
    fastaFiles = ['../../DB/uniprot/uniprot_trembl.fasta.gz']

x = []
flag = 0
dic = {}
#for line in gzip.open('../../DB/uniprot/uniprot_sprot_human.dat.gz', 'rt'):
for file in files:
    for line in gzip.open(file, 'rt'):
        #print (line)
        if line[:2] == 'ID':
            id = line.split()[1]
            dic[id] = spc(id)
            if 'trembl' in file:
                dic[id].uniprotkb = 'TrEMBL'
            else:
                dic[id].uniprotkb = 'Swiss-Prot'
            flag = 0
        elif line[:2] == 'AC' and dic[id].acc == '-':
            dic[id].acc = line.split('AC')[1].replace(' ', '').replace('\n', '').split(';')[0]
        elif line[:2] == 'GN' and dic[id].gene == '-':
            #print (line)
            dic[id].gene = line.split('=')[1].split(' ')[0].replace(' ', '').replace(',', '').replace(';', '').replace('\n', '')
        elif line[:2] == 'KW' and dic[id].kw == 'NO':
            if 'Transmembrane' in line:
                dic[id].kw = 'YES'
        elif line[:2] == 'FT':
            if line.split()[1] in ['SIGNAL']:
                dic[id].signal = 1
            if line.split()[1] in ['SIGNAL', 'TRANSMEM', 'TOPO_DOM']:
                topo = line.split()[1]
                start = line.split()[2].split('.')[0]
                end = line.split()[2].split('.')[-1].replace('\n', '')
                #print (id, line, end.isnumeric())
                if '?' not in start:
                    if int(start) <= 70 and line.split()[1] != 'SIGNAL':
                        flag = 1
                        dic[id].topo += topo+':'+start+':'+end+':'
                    elif line.split()[1] == 'SIGNAL':
                        dic[id].topo += topo+':'+start+':'+end + ': ; '
                        flag = 0
                    else:
                        flag = 0
            else:
                if '/note=' in line.split()[1] and flag == 1:
                    dic[id].topo += line.split()[1].split('"')[1] + '; '
                    if dic[id].type == '-':
                        dic[id].type = line.split()[1].split('"')[1].replace(',', '').replace(';', '')
                    flag = 0
            #break

print (len(dic))
#sys.exit()

print ('Reading fasta')
for files in fastaFiles:
    for line in gzip.open(files, 'rt'):
        if line[0] == '>':
            id = line.split('|')[2].split()[0]
            flag = 0
            if id in dic:
                flag = 1
        else:
            if flag == 1:
                dic[id].fasta += line.replace('\n', '')
'''
for line in open('../data/mechismo_input_uniprot_muts_mods_v9.txt', 'r'):
    id = line.split('\t')[1]
    if id in dic:
        mut = line.split('\t')[0].split('/')[1].split()[0]
        if line.split('\t')[3].split()[0] == 'VARIANT' and '(in'  in line:
            #print (line.split('\t')[3].split('(')[1])
            dis = line.split('\t')[3].split('(')[1]
            if len(dis.split()) >1:
                dis = dis.split()[1]
            dis = dis.replace(';','').replace(':', '')
            if dis.isupper():
                #print (line)
                dic[id].mut += 1
'''

map = {}
mapgid = {}
if species == 'human':
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
else:
    for line in gzip.open('/home/gurdeep/projects/DB/uniprot/CAEEL_6239_idmapping.dat.gz', 'rt'):
        if line.split('\t')[1] == 'UniProtKB-ID':
            acc = line.split('\t')[0]
            id = line.split('\t')[2].replace('\n', '')
            map[acc] = id
        elif line.split('\t')[1] == 'Gene_Name':
            if line.split('\t')[0] == acc:
                gene = line.split('\t')[2].replace('\n', '')
                if gene not in mapgid:
                    mapgid[gene] = id

for line in open('../data/topfind_output.tsv', 'rt'):
    acc = line.split(' ')[0].replace(' ', '')
    position = line.split(',')[2].replace(' ', '')
    if acc in map:
        id = map[acc]
        if position not in dic[id].topfind:
            dic[id].topfind.append(position)

#sys.exit()
#for num, files in enumerate(['clinvar.tsv.gz', 'cosmic.tsv.gz']):
for num, files in enumerate(['clinvar_missense.tsv.gz', 'cosmic.tsv.gz']):
    for line in gzip.open('../data/'+files, 'rt'):
        if line[0] != '#':
            gene = line.split('\t')[0]
            if gene in mapgid:
                id = mapgid[gene]
                if id in dic:
                    if num == 0:
                        #dic[id].clinvar_mut = int(line.split('\t')[1].replace('\n', ''))
                        dic[id].clinvar_mut.append(line.split('\t')[1].replace('\n', ''))
                        dic[id].mut.append(line.split('\t')[1].replace('\n', ''))
                        #dic[id].clinvar_mut.append(line.split('\t')[1].split('.')[1])
                        #dic[id].mut.append(line.split('\t')[1].split('.')[1])
                    else:
                        #dic[id].cosmic_mut = int(line.split('\t')[1].replace('\n', ''))
                        num_sampleid = int(line.split('\t')[2].replace('\n', ''))
                        if num_sampleid >= 3:
                            #dic[id].cosmic_mut += 1
                            dic[id].cosmic_mut.append(line.split('\t')[1].split('.')[1])
                            dic[id].mut.append(line.split('\t')[1].split('.')[1])
'''
for id in dic:
    if dic[id].gene in ['SLC38A2', 'STT3B', 'RHBDF2']:
        print (id, dic[id].gene, dic[id].clinvar_mut, dic[id].cosmic_mut)
sys.exit()
'''
for line in gzip.open('../data/uniprot_2020-05_features_Hsa.tsv.gz', 'rt'):
    if line.split('\t')[1] == 'VARIANT' and '->' in line.split('\t')[5]:
        acc = line.split('\t')[0]
        #id = line.split('\t')[1]
        if acc in map:
            id = map[acc]
            if id in dic:
                #mut = line.split('\t')[0].split('/')[1].split()[0]
                print (line)
                mutation = line.split('\t')[5].split('-')[0]
                mutation += str(line.split('\t')[4])
                mutation += line.split('\t')[5].split('>')[1]
                if 'in'  in line:
                    #dic[id].mut += 1
                    dic[id].uniprot_mut.append(mutation)
                    dic[id].mut.append(mutation)

if species == 'human':
    SIGNALp = '../data/signalp/'
else:
    SIGNALp = '../data/signalp_caeel/'


for id in dic:
    if dic[id].topo == '':
        dic[id].topo = 'NA'
    #print (id, dic[id].acc)
    if os.path.isfile(SIGNALp+id+'_long_noTM.out'):
        #print ('found', id)
        for line in open(SIGNALp+id+'_long_noTM.out', 'r'):
            if line.split()[0] == 'max.' and line.split()[1]=='Y':
                dic[id].prob_notm = float(line.replace('\n','').split()[3])
                dic[id].pos_notm = int(line.replace('\n','').split()[2])
            if line.split()[0] == 'D':
                dic[id].pred_notm = line.replace('\n','').split()[-1]
        for line in open(SIGNALp+id+'_long_TM.out', 'r'):
            if line.split()[0] == 'max.' and line.split()[1]=='Y':
                dic[id].prob_tm = float(line.replace('\n','').split()[3])
                dic[id].pos_tm = int(line.replace('\n','').split()[2])
            if line.split()[0] == 'D':
                dic[id].pred_tm = line.replace('\n','').split()[-1]
    else:
        print (SIGNALp+id+'_long_noTM.out')


df = pd.DataFrame()
data = []
cleavedSeq = ''
fastaAll = ''
l = 'ID\tGENE\tACC\tUniProt-KB\tTM protein\tSignalp\tTotal_MUT\tUniProt_MUT\tClinVar_MUT\tCOSMIC_MUT\tPosition\tSignalP_TM\tSignalP_noTM\tN_terminus\tTopology\tTopfind\tTrypsinCleavageSite\tTrypsinPepLength\tTrypsinPeptide\tLegumainCleavageSite\tLegumainPepLength\tLegumainPeptide\tGluCCleavageSite\tGluCPepLength\tGluCPeptide\n'
for id in dic:
    if dic[id].fasta != '':
        fastaAll += '>' + id + '\n' + dic[id].fasta + '\n'
    #if dic[id].pos_notm == dic[id].pos_tm:
    #if dic[id].type in ['Cytoplasmic']:
    #if dic[id].signal in [0] and dic[id].type == 'Cytoplasmic' and dic[id].kw == 'YES':
    #if dic[id].signal in [0] and dic[id].kw == 'YES':
    #if dic[id].signal in [0] and dic[id].type == 'Cytoplasmic':
    if dic[id].signal in [0, 1]:
        row = []
        #print (id, dic[id].gene, dic[id].acc, dic[id].kw, dic[id].mut, dic[id].prob_tm, dic[id].prob_notm, dic[id].type, dic[id].topo)
        l += id+'\t'+ dic[id].gene+'\t'+ dic[id].acc+'\t'+ dic[id].uniprotkb + '\t' + dic[id].kw+'\t'+ str(dic[id].signal) + '\t' + str(len(list(set(dic[id].mut)))) + '\t' + str(len(dic[id].uniprot_mut))+'\t'+ str(len(dic[id].clinvar_mut))+'\t'+ str(len(dic[id].cosmic_mut))+'\t'+ str(dic[id].pos_notm) + '\t' + str(dic[id].prob_tm)+'\t'+ str(dic[id].prob_notm)+'\t'+ dic[id].type+'\t'+ dic[id].topo+'\t'+ ','.join(dic[id].topfind)+'\n'
        row.append(id)
        row.append(dic[id].gene)
        row.append(dic[id].acc)
        row.append(dic[id].uniprotkb)
        row.append(dic[id].kw)
        row.append(dic[id].signal)
        row.append(len(list(set(dic[id].mut))))
        row.append(len(dic[id].uniprot_mut))
        #row.append(dic[id].clinvar_mut)
        row.append(len(dic[id].clinvar_mut))
        row.append(len(dic[id].cosmic_mut))
        row.append(dic[id].pos_notm)
        row.append(dic[id].prob_tm)
        row.append(dic[id].prob_notm)
        row.append(dic[id].type)
        '''
        if dic[id].type in ['Cytoplasmic', 'Extracellular', 'Lumenal']:
            row.append(dic[id].type)
        else:
            row.append('Other')
        '''
        row.append(dic[id].topo)
        row.append(','.join(dic[id].topfind))
        ## Fetch enzyme cleavage sites
        for enzyme, enzymeName in zip(['[R]', '[DN]', '[DE]'], ['Trypsin', 'Legumain', 'GluC']):
            if dic[id].fasta != '':
                print (id, 'found')
                enzymeCleavageSite = fetchEnzymeCleavageSites(dic[id].fasta, dic[id].pos_notm, enzyme)
                if enzymeCleavageSite == None:
                    row.append('-')
                    row.append('-')
                    row.append('-')
                else:
                    enzymeCleavageSite == int(enzymeCleavageSite)
                    row.append(enzymeCleavageSite)
                    row.append(enzymeCleavageSite - dic[id].pos_notm)
                    row.append(dic[id].fasta[dic[id].pos_notm-1:enzymeCleavageSite-1])
                    cleavedSeq += '>' + id + "|" + enzymeName + '\n' + dic[id].fasta[dic[id].pos_notm-1:enzymeCleavageSite-1] + '\n'
            else:
                print (id)
                row.append('-')
                row.append('-')
                row.append('-')

        '''
        if dic[id].prob_notm >= 0.5 and dic[id].prob_tm <= 0.55:
            row.append('YES')
        else:
            row.append('NO')
        '''
        data.append(row)

#print (cleavedSeq)
open('part_a_cleavedSeq.fasta', 'w').write(cleavedSeq)
open('fastaAll.fasta', 'w').write(fastaAll)
# sys.exit()
df = pd.DataFrame(data, columns=['ID', 'GENE', 'ACC', 'UniProtKB', 'TM protein', 'Signalp', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'SignalP_TM', 'SignalP_noTM', 'N_terminus', 'Topology', 'Topfind', 'TrypsinCleavageSite', 'TrypsinPepLength', 'TrypsinPeptide', 'LegumainCleavageSite','LegumainPepLength', 'LegumainPeptide', 'GluCCleavageSite', 'GluCPepLength', 'GluCPeptide'])
print (df.shape)
'''
fig=pylab.figure()
ax = fig.add_axes([0.1,0.1,0.85,0.85])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES', 'NO'])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES', 'NO'])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES'])
ax = sns.scatterplot(data=df, x="SignalP_noTM", y="SignalP_TM", hue="N_terminus", style='TM protein', legend='brief', style_order=['YES'])
#ax.set_title('Siginificant features', fontsize=7)
ax.set_xticks(np.arange(0,1.1,0.1))
ax.set_yticks(np.arange(0,1.1,0.1))
ax.set_xlabel('SignalP_noTM', fontsize=7)
ax.set_ylabel('SignalP_TM', fontsize=7)
ax.tick_params(axis='y', labelsize=6, left=True, labelleft=True, right=False, labelright=False)
ax.tick_params(axis='x', labelsize=6, bottom=True, labelbottom=True, top=False, labeltop=False)
ax.legend(fontsize=5.5)
ax.grid(True, linewidth=0.001, color='grey', linestyle='--')
ax.plot([0.5,0.5], [0.0,1.0], color='black', linewidth=0.5, linestyle='--')
ax.plot([0.0,1.0], [0.55,0.55], color='black', linewidth=0.5, linestyle='--')
ax.plot([0.0,1.0], [0.0,1.0], color='black', linewidth=0.5, linestyle='--')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
'''
if species != 'human':
    df.to_csv('../data/'+species+'/part_a_data.tsv', sep='\t', index=False)
else:
    df.to_csv('../data/part_a_data.tsv', sep='\t', index=False)
#df.to_csv('../data/try.tsv', sep='\t', index=False)

#plt.savefig('../plots/all.svg')
#plt.savefig('../plots/figA.svg')
#open('../data/part_a_data.tsv', 'w').write(l)
#open('../data/data_try.tsv', 'w').write(l)
#plt.savefig('../plots/all_cytoplasmic.svg')
#plt.savefig('../plots/all_TM.svg')
#plt.savefig('../plots/all_cytoplasmic_TM.svg')
#plt.show()
print ('completed')
