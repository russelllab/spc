import os, sys, requests
import os.path

for line in open('../../data/top_find_candidates.tsv', 'r'):
    acc = line.replace('\n', '')
    if os.path.exists(acc) == False:
        os.system('wget https://topfind.clip.msl.ubc.ca/proteins/show/'+acc+' &')
        '''
        r = requests.get('https://topfind.clip.msl.ubc.ca/proteins/show/'+acc)
        print (acc)
        open(acc, 'w').write(r.text)
        '''
