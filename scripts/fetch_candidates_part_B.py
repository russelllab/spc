## To consider TM_start-5 + x AAs such that the total is 70

## Script to extract candidates from UniProt/SwissProt's human.dat.gz file
import os, sys, gzip, requests

class spc:
    def __init__(self, id):
        self.id = ''
        self.acc = '-'
        self.gene = '-'
        self.kw = 'noTM'
        self.signal = 0
        self.topo = ''
        self.fasta = ''
        self.positions = {}
        self.transmem = {}

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
    elif line[:2] == 'AC' and dic[id].acc == '-':
        dic[id].acc = line.split('AC')[1].replace(' ', '').replace('\n', '').split(';')[0]
    elif line[:2] == 'GN' and dic[id].gene == '-':
        #print (line)
        dic[id].gene = line.split('=')[1].split(' ')[0].replace(' ', '').replace(',', '').replace(';', '').replace('\n', '')
    elif line[:2] == 'KW' and dic[id].kw == 'noTM':
        if 'Transmembrane' in line:
            dic[id].kw = 'TM'
    elif line[:2] == 'FT':
        if line.split()[1] in ['SIGNAL']:
            dic[id].signal = 1
        if line.split()[1] in ['SIGNAL', 'TRANSMEM', 'TOPO_DOM']:
            topo = line.split()[1]
            if topo == 'TOPO_DOM':
                flag = 1
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
                    dic[id].transmem[topo+str(count)] = (type,start,end)
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

print (dic['SIA8F_HUMAN'].transmem)

for line in gzip.open('../../DB/uniprot/uniprot_sprot.fasta.gz', 'rt'):
    if line[0] == '>':
        id = line.split('|')[2].split()[0]
        if id in dic:
            flag = 1
        else:
            flag = 0
    elif flag == 1:
        dic[id].fasta += line.replace('\n', '')

for id in dic:
    if dic[id].fasta == '':
        r = requests.get('https://www.uniprot.org/uniprot/'+dic[id].acc+'.fasta')
        dic[id].fasta = r.text
        print ('Fetched', id)

## Extract Type II TMDs (N-term: Cytoplasmic; C-term: Extracellular)
l = ''
for id in dic:
    if dic[id].fasta != '' and len(dic[id].fasta) >= 70:
        if dic[id].kw == 'TM' and dic[id].signal == 0:
            for tm in dic[id].transmem:
                type, start, end = dic[id].transmem[tm]
                if type == 'Cytoplasmic':
                    name = id+'_'+str(tm)
                    if start <= 5:
                        start_seq = 0
                        end_seq = 70
                    else:
                        start_seq = start-5-1
                        end_seq = start_seq+70
                        if end_seq > len(dic[id].fasta):
                            start_seq = len(dic[id].fasta) - 70
                            end_seq = len(dic[id].fasta)
                    l += '>'+name+'\n'+dic[id].fasta[start_seq:end_seq]+'\n'
                    if len(dic[id].fasta[start_seq:end_seq]) != 70:
                        print (id, name, start_seq, end_seq, len(dic[id].fasta))
                        sys.exit()

#print (l)
#sys.exit()
'''
for id in dic:
    if dic[id].kw == 'noTM' and dic[id].signal == 1:
        #if id == 'HPTR_HUMAN':
        print (id, dic[id].kw, dic[id].topo, dic[id].fasta)
'''
l = ''
z = '#Name\tPEP_START\tPEP_END\tDOM_START\tDOM_END\n'
for id in dic:
    if dic[id].fasta != ''  and len(dic[id].fasta) >= 70:
        for tm in dic[id].transmem:
            type, start, end = dic[id].transmem[tm]
            name = id+'_'+str(tm)
            if start <= 5:
                start_seq = 0
                end_seq = 70
            else:
                start_seq = start-5-1
                end_seq = start_seq+70
                if end_seq > len(dic[id].fasta):
                    start_seq = len(dic[id].fasta) - 70
                    end_seq = len(dic[id].fasta)
            m = '>'+name+'\n'+dic[id].fasta[start_seq:end_seq]
            l += m + '\n'
            #print (name, start-5-1, start_seq+70, start, end)
            z += name + '\t' + str(start-5) + '\t' + str(start_seq+70) + '\t' + str(start) + '\t' +  str(end) + '\n'
            if os.path.isfile('../data/part_b_fasta/'+name+'.fasta') == False:
                open('../data/part_b_fasta/'+name+'.fasta', 'w').write(m)
            if os.path.isfile('../data/part_B_signalp/'+name+'_long_TM.out') == False:
                os.system('/home/gurdeep/Downloads/signalp-4.1/signalp_TM -f long -s best ../data/part_b_fasta/'+name+'.fasta > ../data/part_B_signalp/'+name+'_long_TM.out')
            if os.path.isfile('../data/part_B_signalp/'+name+'_long_noTM.out') == False:
                os.system('/home/gurdeep/Downloads/signalp-4.1/signalp_noTM -f long -s best ../data/part_b_fasta/'+name+'.fasta > ../data/part_B_signalp/'+name+'_long_noTM.out')
            #sys.exit()

#open('../data/uniprot_sprot_human_signalp_input_part_b.fasta', 'w').write(l)
open('../data/part_b_annotations.tsv', 'w').write(z)
