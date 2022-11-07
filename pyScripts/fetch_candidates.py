## Script to extract candidates from
## UniProt/SwissProt's human.dat.gz file
import os, gzip, sys

## Select the species
species = 'caeel'
#species = 'human'
sprot = 1
trembl = 1

if sprot == 1 and trembl == 1:
    #files = ['../../DB/uniprot/uniprot_sprot_human.dat.gz', '../../DB/uniprot/uniprot_trembl_human.dat.gz']
    files = ['../../DB/uniprot/uniprot_sprot_caeel.dat.gz', '../../DB/uniprot/uniprot_trembl_caeel.dat.gz']
    fastaFiles = ['../../DB/uniprot/uniprot_sprot.fasta.gz', '../../DB/uniprot/uniprot_trembl.fasta.gz']
elif sprot == 1 and trembl == 0:
    #files = ['../../DB/uniprot/uniprot_sprot_human.dat.gz']
    files = ['../../DB/uniprot/uniprot_sprot_caeel.dat.gz']
    fastaFiles = ['../../DB/uniprot/uniprot_sprot.fasta.gz']
elif sprot == 0 and trembl == 1:
    #files = ['../../DB/uniprot/uniprot_trembl_human.dat.gz']
    files = ['../../DB/uniprot/uniprot_trembl_caeel.dat.gz']
    fastaFiles = ['../../DB/uniprot/uniprot_trembl.fasta.gz']

class spc:
    def __init__(self, id):
        self.id = ''
        self.kw = 'noTM'
        self.signal = 0
        self.topo = ''
        self.fasta = ''

x = []
flag = 0
dic = {}
for file in files:
    for line in gzip.open(file, 'rt'):
    #for line in gzip.open('../../DB/uniprot/uniprot_sprot_human.dat.gz', 'rt'):
        #print (line)
        if line[:2] == 'ID':
            id = line.split()[1]
            dic[id] = spc(id)
            flag = 0
        elif line[:2] == 'KW' and dic[id].kw == 'noTM':
            if 'Transmembrane' in line:
                dic[id].kw = 'TM'
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
                        dic[id].topo += topo+':'+start+':'+end + ':\t'
                        flag = 0
                    else:
                        flag = 0
            else:
                if '/note=' in line.split()[1] and flag == 1:
                    dic[id].topo += line.split()[1].split('"')[1] + '\t'
                    flag = 0
            #break

print (len(dic))

for file in fastaFiles:
    for line in gzip.open(file, 'rt'):
        if line[0] == '>':
            id = line.split('|')[2].split()[0]
            if id in dic:
                flag = 1
            else:
                flag = 0
        elif flag == 1:
            if len(dic[id].fasta) < 70:
                for char in line.replace('\n', ''):
                    dic[id].fasta += char
                    #if id == 'YO002_HUMAN':
                    #    print (dic[id].fasta, len(dic[id].fasta))
                    if len(dic[id].fasta) == 70:
                        flag = 0
                        break

print (len(dic))
#sys.exit()
'''
for id in dic:
    if dic[id].kw == 'noTM' and dic[id].signal == 1:
        #if id == 'HPTR_HUMAN':
        print (id, dic[id].kw, dic[id].topo, dic[id].fasta)
'''
l = ''
for id in dic:
    if dic[id].fasta != '':
        m = '>'+id + '\n' + dic[id].fasta + '\n'
        l += m
        open('../data/fasta_'+species+'/'+id+'.fasta', 'w').write(m)
        os.system('/home/gurdeep/Downloads/signalp-4.1/signalp_TM -f long -s best ../data/fasta_'+species+'/'+id+'.fasta > ../data/signalp_'+species+'/'+id+'_long_TM.out')
        os.system('/home/gurdeep/Downloads/signalp-4.1/signalp_noTM -f long -s best ../data/fasta_'+species+'/'+id+'.fasta > ../data/signalp_'+species+'/'+id+'_long_noTM.out')

#open('../data/uniprot_sprot_human_signalp_input.fasta', 'w').write(l)
open('../data/uniprot_'+species+'_signalp_input.fasta', 'w').write(l)
