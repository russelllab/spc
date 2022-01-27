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
        self.subcellular = []
        self.subcell = 'Unannotated in UniProt'
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
        self.mut = 0

x = []
flag = 0
dic = {}
for line in gzip.open('../../DB/uniprot/uniprot_sprot_human.dat.gz', 'rt'):
    #print (line)
    if line[:2] == 'ID':
        id = line.split()[1]
        dic[id] = spc(id)
        flag = 0
        subcell_flag = 0
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
    elif line[:2] == 'CC':
        if 'SUBCELLULAR LOCATION' in line and line.split()[1] == '-!-':
            #print (id, line.split())
            dic[id].subcellular.append(line.replace('\n', '').lower())
            subcell_flag = 1
            #sys.exit()
        elif line.split()[1] == '-!-':
            subcell_flag = 0
            #sys.exit()
        elif subcell_flag == 1:
            #print (id, line.split())
            dic[id].subcellular.append(line.replace('\n', '').lower())

for id in dic:
    for location in dic[id].subcellular:
        if 'peripheral membrane' in location:
            dic[id].subcell = 'Peripheral membrane proteins'
            break
        elif 'mitochondrion' in location:
            dic[id].subcell = 'Mitochondrion proteins'
            break
        elif 'cytoplasm' in location or 'cytosol' in location or 'nuclear' in location or 'nucleus' in location or 'endoplasmic reticulum' in location:
            dic[id].subcell = 'Soluble (Cytoplasm/Nucleus)'
            break
    if dic[id].subcell == '':
        dic[id].subcell = 'Others'

for id in dic:
    if dic[id].subcell == 'Peripheral membrane proteins' and dic[id].kw == 'NO':
        print (id, dic[id].subcellular)

#sys.exit()
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
l = 'ID\tGENE\tACC\tTM protein\tSignal peptide\tNum MUT\tSub-cellular location\tPosition\tSignalP_TM\tSignalP_noTM\tN-terminus\tTopology\n'
for id in dic:
    #if dic[id].pos_notm == dic[id].pos_tm:
    #if dic[id].type in ['Cytoplasmic']:
    #if dic[id].signal in [0] and dic[id].type == 'Cytoplasmic' and dic[id].kw == 'YES':
    if dic[id].signal in [0] and dic[id].kw == 'YES':
    #if dic[id].signal in [0] and dic[id].type == 'Cytoplasmic':
    #if dic[id].signal in [0] and dic[id].type in ['Cytoplasmic', 'Extracellular', 'Lumenal']:
    #if dic[id].signal in [0]:
        row = []
        dic_sp = {1:'YES', 0:'NO'}
        #print (id, dic[id].gene, dic[id].acc, dic[id].kw, dic[id].mut, dic[id].prob_tm, dic[id].prob_notm, dic[id].type, dic[id].topo)
        l += id+'\t'+dic[id].gene+'\t'+dic[id].acc+'\t'+dic[id].kw+'\t'+str(dic_sp[dic[id].signal])+'\t'+str(dic[id].mut)+'\t'+str(dic[id].subcell)+'\t'+str(dic[id].pos_notm)+'\t'+str(dic[id].prob_tm)+'\t'+str(dic[id].prob_notm)+'\t'+dic[id].type+'\t'+dic[id].topo+'\n'
        row.append(id)
        row.append(dic[id].kw)
        row.append(dic[id].signal)
        row.append(dic[id].mut)
        if dic[id].type in ['Cytoplasmic', 'Extracellular', 'Lumenal']:
            row.append(dic[id].type)
        else:
            row.append('Undefined')
        row.append(dic[id].prob_tm)
        row.append(dic[id].prob_notm)
        '''
        if dic[id].prob_notm >= 0.5 and dic[id].prob_tm <= 0.55:
            row.append('YES')
        else:
            row.append('NO')
        '''
        if dic[id].signal == 1:
            row.append('Canonical Signal Peptide')
        elif dic[id].kw == 'YES':
            if dic[id].type in ['Cytoplasmic']:
                row.append('Type II (N-t: '+dic[id].type+')')
            elif dic[id].type in ['Extracellular', 'Lumenal']:
                row.append('Type I (N-t: '+'Extracellular/Lumenal)')
            else:
                row.append('Type Undefined')
                #print (id)
                #sys.exit()
        else:
            row.append(dic[id].subcell)

        data.append(row)

#print (data)
df = pd.DataFrame(data, columns=['ID', 'TM protein', 'SIGNALp', 'MUT', 'N-terminus', 'signalp_TM', 'signalp_noTM', 'Categories'])
ho = np.sort(list(set(df['Categories'].to_numpy())))
print (ho)
#sys.exit()
print (df.shape)
sizes = [4 for i in range(0, len(df['Categories'].to_numpy()))]
fig=pylab.figure()
ax = fig.add_axes([0.1,0.1,0.85,0.85])
#colors = ['lightblue', 'pink', 'orange', 'magenta', 'red', 'green', 'blue', 'brown']
#colors = ['pink', 'orange', 'magenta', 'red', 'green', 'blue', 'brown']
colors = ['red', 'green', 'blue']
ax = sns.scatterplot(data=df, palette=colors, x="signalp_noTM", y="signalp_TM", hue="Categories", hue_order=ho, style='TM protein', legend='brief', style_order=['YES', 'NO'], size='MUT')
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES', 'NO'])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES'])
#ax = sns.scatterplot(data=df, x="signalp_noTM", y="signalp_TM", hue="N-terminus", style='TM protein', legend='brief', style_order=['YES'])
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
#open('../data/data_customized.tsv', 'w').write(l)
#plt.show()
#sys.exit()
#plt.savefig('../plots/customized_ppt.svg')
#plt.savefig('../plots/customized_ppt.eps')
#plt.savefig('../plots/customized_ppt.png')
#plt.savefig('../plots/customized_without_signal_peptide_ppt.svg')
#plt.savefig('../plots/customized_without_signal_peptide_ppt.eps')
#plt.savefig('../plots/customized_without_signal_peptide_ppt.png')
#plt.savefig('../plots/customized_without_signal_peptide_mit_ppt.svg')
#plt.savefig('../plots/customized_without_signal_peptide_mit_ppt.eps')
#plt.savefig('../plots/customized_without_signal_peptide_mit_ppt.png')
plt.savefig('../plots/customized_only_TM.svg')
plt.savefig('../plots/customized_only_TM.eps')
plt.savefig('../plots/customized_only_TM.png')
