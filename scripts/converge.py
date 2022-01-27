import os, sys, gzip

class proteins:
    def __init__(self, id):
        self.id = id
        self.topo = {}
        self.mut = []
        self.uniprot_mut = []
        self.cosmic_mut = []
        self.clinvar_mut = []

dic = {}
for line in open('../data/loop_length_mutations.tsv', 'r'):
    if line.split('\t')[1] != 'ID':
        id = line.split('\t')[1]
        topo_dom = line.split('\t')[5]
        topo_dom_start = line.split('\t')[6]
        topo_dom_end = line.split('\t')[7]
        if id not in dic:
            dic[id] = proteins(id)
        dic[id].topo[topo_dom] = (topo_dom_start, topo_dom_end)

print (dic['VTM2B_HUMAN'].topo)

################ Mutation data #################################################
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
################################################################################

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

l = ''
for line in open('../data/loop_mutaions.tsv', 'r'):
    if line.split('\t')[0] != 'ID':
        id = line.split('\t')[0]
        if id in dic:
            precede_topo_dom = line.split('\t')[16]
            if precede_topo_dom != '':
                num = int(precede_topo_dom.split('TOPO_DOM')[1])
                succeed_topo_dom = 'TOPO_DOM' + str(num+1)
                if succeed_topo_dom in dic[id].topo:
                    start, end = dic[id].topo[succeed_topo_dom]
                    l += line.replace('\n', '\t') + succeed_topo_dom + '\t' + str(start) + '\t' + str(end) + '\t'
                    l += str(count_mutations(id, int(start), int(end), dic[id].mut)) + '\t'
                    l += str(count_mutations(id, int(start), int(end), dic[id].uniprot_mut)) + '\t'
                    l += str(count_mutations(id, int(start), int(end), dic[id].clinvar_mut)) + '\t'
                    l += str(count_mutations(id, int(start), int(end), dic[id].cosmic_mut)) + '\n'
                    #print (line, succeed_topo_dom, num)
    else:
        l += line.replace('\n', '\t') + 'SUC_TOP_DOM\tSUC_START\tSUC_END\tSUC_MUT\tSUC_UniProt\tSUC_ClinVar\tSUC_COSMIC\n'

open('../data/converge.tsv', 'w').write(l)
