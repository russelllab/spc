import xmltojson
import json
import requests, os

for files in os.listdir('.'):
    if '.' not in files:
        #print (files)
        for line in open(files, 'r'):
            if 'var domainElem' in line and 'Cleavages' in line:
                #print (line.replace('],[', '\n').replace('[','').replace(']',''))
                for text in line.split('],['):
                    if 'Cleavages' in text:
                        print (files, text)
