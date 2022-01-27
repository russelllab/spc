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
        self.kw = 'NO'
        self.signal = 0
        self.topo = ''
        self.fasta = ''
        self.type = 'NA'
        self.positions = {}
        self.transmem = {}
        self.pos_tm = 0
        self.pos_notm = 0
        self.prob_tm = 0
        self.prob_notm = 0
        self.pred_tm = ''
        self.pred_notm = ''
        self.mut = []
        self.uniprot_mut = []
        self.cosmic_mut = []
        #self.clinvar_mut = 0
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
        count_topo_dom = 0
        topo_dom = None
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
    elif line[:2] == 'KW' and dic[id].kw == 'NO':
        if 'Transmembrane' in line:
            dic[id].kw = 'YES'
    elif line[:2] == 'FT':
        if line.split()[1] in ['SIGNAL']:
            dic[id].signal = 1
        if line.split()[1] in ['SIGNAL', 'TRANSMEM', 'TOPO_DOM']:
            topo = line.split()[1]
            if topo == 'TOPO_DOM':
                flag = 1
                count_topo_dom += 1
                topo_dom = topo + str(count_topo_dom)
                continue
            elif topo == 'TRANSMEM':
                count += 1
                start = line.split()[2].split('.')[0]
                end = line.split()[2].split('.')[-1].replace('\n', '')
                if '?' not in start and '?' not in end:
                    start = int(start)
                    end = int(end)
                    flag = 0
                    #if id == 'CXB1_HUMAN':
                    #print (count, id, line, type)
                    dic[id].transmem[topo+str(count)] = (topo_dom, type,start,end)
            '''
            if '?' not in start:
                if int(start) <= 70 and line.split()[1] != 'SIGNAL':
                    flag = 1
                    dic[id].topo += topo+':'+start+':'+end+':'
                elif line.split()[1] == 'SIGNAL':
                    dic[id].topo += topo+':'+start+':'+end + ':\t'
                    flag = 0
                else:
                    flag = 0
                '''
        else:
            '''
            if '/note=' in line.split()[1] and flag == 1:
                dic[id].topo += line.split()[1].split('"')[1] + '\t'
                flag = 0
            '''
            if '/note=' in line.split()[1] and flag == 1:
                type = line.split()[1].split('"')[1]
                flag = 0

        #break
    if sq == 1:
        dic[id].fasta += line.replace('\n', '').replace(' ', '')

#print (dic['PMP22_HUMAN'].transmem)
#sys.exit()
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

#print (dic['CXB1_HUMAN'].mut)
#sys.exit()
## Store START of peptide
map_start = {}
for line in open('../data/part_c_annotations.tsv', 'r'):
    if line[0] != '#':
        name = line.split('\t')[0]
        start = int(line.split('\t')[1])
        map_start[name] = start

for line in open('../data/topfind_output.tsv', 'rt'):
    acc = line.split(' ')[0].replace(' ', '')
    position = line.split(',')[2].replace(' ', '')
    id = map[acc]
    if position not in dic[id].topfind:
        dic[id].topfind.append(position)


df = pd.DataFrame()
data = []
#l = 'ID\tGENE\tACC\tTM protein\tNum MUT\tPosition\tSignalP_TM\tSignalP_noTM\tN-terminus\tTopology\n'
l = 'ID\tGENE\tACC\tTM protein\tTotal_MUT\tUniProt MUT\tClinVar MUT\tCOSMIC MUT\tPosition\tSignalP_TM\tSignalP_noTM\tN-terminus\tTopology\tTopfind\n'
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
        if dic[id].topo == '':
            dic[id].topo = 'NA'
        #print (id, dic[id].acc)
        for tm in dic[id].transmem:
            total_count += 1
            topo_dom, type, start, end = dic[id].transmem[tm]
            if type == 'Cytoplasmic':
                count_type2 += 1
            elif type in ['Extracellular', 'Lumenal']:
                count_type1 += 1
                #continue
            else:
                count_type3 += 1
                #print (type)
                #continue

            if start <= 70 and int(tm.split('TRANSMEM')[1]) == 1:
                #print (id, tm, dic[id].transmem)
                #sys.exit()
                count_dom1 += 1
                continue
            name = id+'_'+str(tm)
            #if 'PMP22_HUMAN' in name:
            #    print (name, type, start, end)
            if True:
                start_seq = start-1
                end_seq = start_seq+70
                if end_seq > len(dic[id].fasta):
                    #start_seq = len(dic[id].fasta) - 70
                    end_seq = len(dic[id].fasta)
                    '''
                    check = []
                    for tm_check in dic[id].transmem:
                        type_check, start_check, end_check = dic[id].transmem[tm_check]
                        if 'PMP22_HUMAN' in name:
                            #print (dic[id].transmem)
                            print (type_check, start_check, end_check, start_seq, end_seq, len(dic[id].fasta))
                        if (start_check<=end_seq and start_check>=start_seq) or (end_check<=end_seq and end_check>=start_seq):
                            row =[]
                            #print (type_check, start_check, end_check, end='\t')
                            row.append(type_check)
                            row.append(start_check)
                            row.append(end_check)
                            check.append(row)
                    #if len(check) == 1 and type != 'Unknown':
                    if 'PMP22_HUMAN' in name:
                        print (check)
                    if len(check) > 1:
                        count_excl += 1
                        if 'PMP22_HUMAN' in name:
                            print (name)
                        #print (id, check)
                        continue
                        sys.exit()
                    '''

            prob_tm = 0.0
            prob_notm = 0.0
            if os.path.isfile('../data/part_c_signalp/'+name+'_long_noTM.out'):
                #print ('YES')
                for line in open('../data/part_c_signalp/'+name+'_long_noTM.out', 'r'):
                    if line.split()[0] == 'max.' and line.split()[1]=='Y':
                        prob_notm = float(line.replace('\n','').split()[3])
                        pos_notm = int(line.replace('\n','').split()[2]) + map_start[name] - 1
                        dic[id].pos_notm = pos_notm
                    if line.split()[0] == 'D':
                        pred_notm = line.replace('\n','').split()[-1]
                for line in open('../data/part_c_signalp/'+name+'_long_TM.out', 'r'):
                    if line.split()[0] == 'max.' and line.split()[1]=='Y':
                        prob_tm = float(line.replace('\n','').split()[3])
                        pos_tm = int(line.replace('\n','').split()[2]) + map_start[name] - 1
                        dic[id].pos_tm = pos_tm
                    if line.split()[0] == 'D':
                        pred_tm = line.replace('\n','').split()[-1]
            else:
                print ('No', '../data/part_c_signalp/'+name+'_long_noTM.out')
                #sys.exit()

            row = []
            #print (name, dic[id].gene, dic[id].acc, dic[id].kw, dic[id].mut, prob_tm, prob_notm, dic[id].type, dic[id].topo)
            #l += id+'\t'+ dic[id].gene+'\t'+ dic[id].acc+'\t'+ dic[id].kw+'\t'+ str(dic[id].mut)+'\t'+ str(dic[id].pos_notm) + '\t' + str(dic[id].prob_tm)+'\t'+ str(dic[id].prob_notm)+'\t'+ dic[id].type+'\t'+ dic[id].topo+'\n'
            row.append(id)
            row.append(dic[id].gene)
            row.append(dic[id].acc)
            row.append(name)
            row.append(tm)
            row.append(start)
            row.append(end)
            row.append(dic[id].pos_notm)
            row.append(dic[id].kw)
            row.append(dic[id].signal)
            row.append(len(list(set(dic[id].mut))))
            row.append(len(dic[id].uniprot_mut))
            #row.append(dic[id].clinvar_mut)
            row.append(len(dic[id].clinvar_mut))
            #row.append(len(list(set(dic[id].clinvar_mut))))
            row.append(len(dic[id].cosmic_mut))
            row.append(5)
            #row.append(type)

            if type in ['Cytoplasmic', 'Extracellular', 'Lumenal', 'Nuclear']:
                row.append(type)
            else:
                row.append('Other')

            row.append(topo_dom)

            row.append(prob_tm)
            row.append(prob_notm)
            row.append(','.join(dic[id].topfind))
            if prob_notm >= 0.5 and prob_tm <= 0.5 and type == 'Cytoplasmic':
                row.append('YES')
            else:
                row.append('NO')
            #if 'PMP22_HUMAN' in name:
            #    print (row)
            data.append(row)

#sys.exit()
df = pd.DataFrame(data, columns=['ID', 'GENE', 'ACC', 'Name', 'TMD', 'Start', 'End', 'Position', 'TM protein', 'Signalp', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Size', 'N_terminus', 'TOPO_DOM', 'SignalP_TM', 'SignalP_noTM', 'Topfind', 'Highlight'])
print (df.shape, count_type1, count_type2, count_type3, total_count, count_dom1, count_excl)
fig=pylab.figure()
ax = fig.add_axes([0.1,0.1,0.85,0.85])
ax = sns.scatterplot(data=df, x="SignalP_noTM", y="SignalP_TM", hue="N_terminus", style='Signalp', legend='brief', s=15)
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
#plt.savefig('../part_C_plots/fig_partc.eps')
#plt.savefig('../part_C_plots/fig_partc.svg')
#plt.savefig('../part_C_plots/fig_partc.png')
df.to_csv('../data/loop_mutaions.tsv', sep='\t', index=False)
#df.to_csv('../data/part_c_data_try.tsv', sep='\t'. index=False)
#plt.savefig('../part_C_plots/all_cytoplasmic.svg')
#plt.savefig('../part_C_plots/all_TM.svg')
#plt.savefig('../part_C_plots/all_cytoplasmic_TM.svg')
#plt.show()
print ('completed')
