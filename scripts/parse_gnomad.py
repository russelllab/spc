import os, sys, gzip

class ensembl:
    def __init__(self, transcript, acc):
        self.transcript = transcript
        self.acc = acc
        self.mutations = []

dic = {}
for line in gzip.open('../../DB/uniprot/HUMAN_9606_idmapping.dat.gz', 'rt'):
    if line.split('\t')[1] == 'Ensembl_TRS':
        transcript = line.split('\t')[2].replace('\n', '')
        acc = line.split('\t')[0].replace('\n', '')
        if transcript not in dic:
            dic[transcript] = ensembl(transcript, acc)

dir = '../data/'
count = 0
for files in os.listdir(dir):
    if files.startswith('gnomad'):
        print (files)
        for line in gzip.open(dir+'gnomad.exomes.r2.1.1.sites.Y_parsed.txt.gz', 'rt'):
            if 'UniProt_ID' not in line:
                acc = line.split('\t')[0].replace(' ', '')
                gene = line.split('\t')[1].replace(' ', '')
                transcript = line.split('\t')[3].split('.')[0].replace(' ', '')
                mutation = line.split('\t')[4].replace(' ', '')
                variant = line.split('\t')[5].replace(' ', '')
                filter = line.split('\t')[6].replace(' ', '')
                if variant.split()[0] == 'missense_variant' and filter.split()[0] == 'PASS':
                    if transcript in dic:
                        dic[transcript].mutations.append(mutation)
                        #if count == 100:
                        #    break

l = ''
for transcript in dic:
    if dic[transcript].mutations != []:
        #print (dic[transcript].mutations)
        l += dic[transcript].acc + '\t' + str(transcript) + '\t' + str(len(list(set(dic[transcript].mutations)))) + '\n'
print (l)
