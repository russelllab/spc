## Script to plot SPC
import os, sys, gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pylab
import numpy as np

TMCutoff = 0.6
noTMCutoff = 0.5

class spc:
    def __init__(self, id):
        self.id = ''
        self.acc = '-'
        self.gene = '-'
        self.kw = 'NO'
        self.signal = 0
        self.topo = ''
        self.fasta = ''
        self.type = 'NA'
        self.pos_tm = 0
        self.pos_notm = 0
        self.prob_tm = 0
        self.prob_notm = 0
        self.pred_tm = ''
        self.pred_notm = ''
        self.mut = []
        self.uniprot_mut = []
        self.cosmic_mut = []
        #self.clinvar_mut = []
        self.clinvar_mut = 0
        self.topfind = []

x = []
flag = 0
dic = {}
for line in gzip.open('../../DB/uniprot/uniprot_sprot_human.dat.gz', 'rt'):
    #print (line)
    if line[:2] == 'ID':
        id = line.split()[1]
        dic[id] = spc(id)
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
                if dic[id].type == 'NA':
                    dic[id].type = line.split()[1].split('"')[1].replace(',', '').replace(';', '')
                flag = 0
        #break
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

for line in open('../data/topfind_output.tsv', 'rt'):
    acc = line.split(' ')[0].replace(' ', '')
    position = line.split(',')[2].replace(' ', '')
    id = map[acc]
    if position not in dic[id].topfind:
        dic[id].topfind.append(position)
#sys.exit()
for num, files in enumerate(['clinvar.tsv.gz', 'cosmic.tsv.gz']):
    for line in gzip.open('../data/'+files, 'rt'):
        if line[0] != '#':
            gene = line.split('\t')[0]
            if gene in mapgid:
                id = mapgid[gene]
                if id in dic:
                    if num == 0:
                        dic[id].clinvar_mut = int(line.split('\t')[1].replace('\n', ''))
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
        id = map[acc]
        if id in dic:
            #mut = line.split('\t')[0].split('/')[1].split()[0]
            #print (line)
            mutation = line.split('\t')[5].split('-')[0]
            mutation += str(line.split('\t')[4])
            mutation += line.split('\t')[5].split('>')[1]
            if 'in'  in line:
                #dic[id].mut += 1
                dic[id].uniprot_mut.append(mutation)
                dic[id].mut.append(mutation)

for id in dic:
    if dic[id].topo == '':
        dic[id].topo = 'NA'
    #print (id, dic[id].acc)
    if os.path.isfile('../data/signalp/'+id+'_noTM.out'):
        for line in open('../data/signalp/'+id+'_noTM.out', 'r'):
            if line.split()[0] == 'max.' and line.split()[1]=='Y':
                dic[id].prob_notm = float(line.replace('\n','').split()[3])
                dic[id].pos_notm = int(line.replace('\n','').split()[2])
            if line.split()[0] == 'D':
                dic[id].pred_notm = line.replace('\n','').split()[-1]
        for line in open('../data/signalp/'+id+'_TM.out', 'r'):
            if line.split()[0] == 'max.' and line.split()[1]=='Y':
                dic[id].prob_tm = float(line.replace('\n','').split()[3])
                dic[id].pos_tm = int(line.replace('\n','').split()[2])
            if line.split()[0] == 'D':
                dic[id].pred_tm = line.replace('\n','').split()[-1]

df = pd.DataFrame()
data = []
l = 'ID\tGENE\tACC\tTM protein\tTotal_MUT\tUniProt MUT\tClinVar MUT\tCOSMIC MUT\tPosition\tSignalP_TM\tSignalP_noTM\tN-terminus\tTopology\tTopfind\n'
for id in dic:
    #if dic[id].pos_notm == dic[id].pos_tm:
    #if dic[id].type in ['Cytoplasmic']:
    #if dic[id].signal in [0] and dic[id].type == 'Cytoplasmic' and dic[id].kw == 'YES':
    #if dic[id].signal in [0] and dic[id].kw == 'YES':
    #if dic[id].signal in [0] and dic[id].type == 'Cytoplasmic':
    if dic[id].signal in [0, 1]:
        row = []
        #print (id, dic[id].gene, dic[id].acc, dic[id].kw, dic[id].mut, dic[id].prob_tm, dic[id].prob_notm, dic[id].type, dic[id].topo)
        l += id+'\t'+ dic[id].gene+'\t'+ dic[id].acc+'\t'+ dic[id].kw+'\t'+ str(len(list(set(dic[id].mut)))) + '\t' + str(len(dic[id].uniprot_mut))+'\t'+ str(dic[id].clinvar_mut)+'\t'+ str(len(dic[id].cosmic_mut)) +'\t'+ str(dic[id].pos_notm) + '\t' + str(dic[id].prob_tm)+'\t'+ str(dic[id].prob_notm)+'\t'+ dic[id].type+'\t'+ dic[id].topo+'\t'+ ','.join(dic[id].topfind)+'\n'
        row.append(id)
        row.append(dic[id].kw)
        row.append(dic[id].signal)
        row.append(len(list(set(dic[id].mut))))
        row.append(len(dic[id].uniprot_mut))
        row.append(dic[id].clinvar_mut)
        #row.append(len(list(set(dic[id].clinvar_mut))))
        row.append(len(dic[id].cosmic_mut))
        if dic[id].type in ['Cytoplasmic', 'Extracellular', 'Lumenal']:
            row.append(dic[id].type)
        else:
            row.append('Other')
        row.append(dic[id].prob_tm)
        row.append(dic[id].prob_notm)
        row.append(','.join(dic[id].topfind))
        if dic[id].gene in ['GJB2', 'GJB1', 'GJB4', 'SYVN1', 'PMP22']:
            row.append(3)
            row.append(5)
        elif dic[id].prob_notm >= noTMCutoff and dic[id].prob_tm <= TMCutoff and dic[id].type == 'Cytoplasmic':
            row.append(2)
            row.append(3)
        elif dic[id].signal in [1]:
            row.append(0)
            row.append(3)
        else:
            row.append(1)
            row.append(3)
        data.append(row)
#sys.exit()
df = pd.DataFrame(data,
                    columns=['ID',
                            'TM protein',
                            'SIGNALp',
                            'Total_MUT',
                            'Uniprot_MUT',
                            'ClinVar_MUT',
                            'COSMIC_MUT',
                            'N-terminus',
                            'signalp_TM',
                            'signalp_noTM',
                            'Topfind',
                            'Highlight',
                            'Size'])
df = df.sort_values(by=['Highlight'], ascending=True)
print (df.shape)
fig=pylab.figure()
ax = fig.add_axes([0.1,0.1,0.85,0.85])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES', 'NO'])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES', 'NO'])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES'])
ax = sns.scatterplot(data=df[df["Highlight"].isin([0,1])],
                    x="signalp_noTM",
                    y="signalp_TM",
                    hue="Highlight",
                    style='SIGNALp',
                    size='Size',
                    #legend='brief',
                    legend=False,
                    hue_order=[0,1,2,3],
                    palette=['mediumturquoise', 'goldenrod', 'mediumorchid', 'purple'],
                    sizes=(75, 150),
                    alpha=0.4)

sns.scatterplot(data=df[df["Highlight"].isin([2,3])],
                    ax = ax,
                    x="signalp_noTM",
                    y="signalp_TM",
                    hue="Highlight",
                    style='SIGNALp',
                    size='Size',
                    #legend='brief',
                    legend=False,
                    hue_order=[0,1,2,3],
                    palette=['mediumturquoise', 'goldenrod', 'mediumorchid', 'purple'],
                    sizes=(75, 150),
                    alpha=1.0)
#ax.set_title('Siginificant features', fontsize=7)
ax.set_xticks(np.arange(0,1.1,0.25))
ax.set_yticks(np.arange(0,1.1,0.25))
ax.set_xlabel('SignalP_noTM', fontsize=12)
ax.set_ylabel('SignalP_TM', fontsize=12)
ax.tick_params(axis='y', labelsize=12, left=True, labelleft=True, right=False, labelright=False)
ax.tick_params(axis='x', labelsize=12, bottom=True, labelbottom=True, top=False, labeltop=False)
#ax.legend(fontsize=5.5)
#ax.grid(True, linewidth=0.001, color='grey', linestyle='--')
ax.plot([noTMCutoff,noTMCutoff], [0.0,1.0], color='purple', linewidth=0.5, linestyle='--')
ax.plot([0.0,1.0], [TMCutoff, TMCutoff], color='purple', linewidth=0.5, linestyle='--')
ax.plot([0.0,1.0], [0.0,1.0], color='black', linewidth=0.5, linestyle='--')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#plt.savefig('../plots/all.svg')
plt.savefig('../plots/figA.svg')
#plt.savefig('../plots/figA.eps')
plt.savefig('../plots/figA.png')
#open('../data/data.tsv', 'w').write(l)
#open('../data/data_try.tsv', 'w').write(l)
#plt.savefig('../plots/all_cytoplasmic.svg')
#plt.savefig('../plots/all_TM.svg')
#plt.savefig('../plots/all_cytoplasmic_TM.svg')
#plt.show()
