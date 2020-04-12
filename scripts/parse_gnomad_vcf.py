# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 21:32:53 2018

@author: anton
"""
from collections import defaultdict
import csv
import re


#nfe - Non-finish European
#fin - Finish European
#eas - East Asian

#nfe_nwe -North-western European
#nfe_bgr - Bulgarian
#nfe_est - Estonian
#nfe_seu - Southern European
#nfe_swe - Swedish
#nfe_onf - other non-Finnish

def find_dat(pref,row):
    try:
        ret_list=[re.compile(r';AN_{}=[a-zA-Z0-9.]*'.format(pref)).findall(row[7])[0].replace(';','').split('=')[1],
            re.compile(r';AC_{}=[a-zA-Z0-9.]*'.format(pref)).findall(row[7])[0].replace(';','').split('=')[1]]
    
        if len(ret_list)==2:
            return(ret_list + [int(ret_list[1])/int(ret_list[0])])
        else:
            return(['0','0','0'])
    except:
        return(['0','0','0'])

#prefs=['nfe','fin','eas','nfe_nwe','nfe_bgr','nfe_est','nfe_seu','nfe_swe','nfe_onf']
prefs=['nfe','fin','eas','afr','sas','nfe_bgr','nfe_swe','nfe_nwe','nfe_est','nfe_seu','nfe_onf']
row_dict=defaultdict(list)
with open('/home/anton/gnomad_ACE2_genomes_passed.vcf','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if ':p' in row[7].split('vep=')[1].split('|')[11]:
            prot = row[7].split('vep=')[1].split('|')[11].split(':')[1]
        else:
            prot = '-'

        ad_row=[row[1],row[2],row[3],row[4],row[7].split('vep=')[1].split('|')[1],row[7].split('vep=')[1].split('|')[2]]
        it=[find_dat(name, row) for name in prefs]
        for frq in it:
            ad_row.extend(frq)
        row_dict[ad_row[0]+ad_row[2]+ad_row[3]]=ad_row

count_inds=[i for i in range(8,39,3)]


with open('/home/anton/gnomad_ACE2_passed_exomes.vcf','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if ':p' in row[7].split('vep=')[1].split('|')[11]:
            prot = row[7].split('vep=')[1].split('|')[11].split(':')[1]
        else:
            prot = '-'

        ad_row=[row[1],row[2],row[3],row[4],row[7].split('vep=')[1].split('|')[1],row[7].split('vep=')[1].split('|')[2]]
        it=[find_dat(name, row) for name in prefs]
        for frq in it:
            ad_row.extend(frq)
        if ad_row[0]+ad_row[2]+ad_row[3] in row_dict.keys():
            for i in range(len(ad_row)):
                if i>5:
                    if i not in count_inds:
                        row_dict[ad_row[0]+ad_row[2]+ad_row[3]][i]=int(row_dict[ad_row[0]+ad_row[2]+ad_row[3]][i]) + int(ad_row[i])
                    else:
                        row_dict[ad_row[0]+ad_row[2]+ad_row[3]][i]=int(row_dict[ad_row[0]+ad_row[2]+ad_row[3]][i-1])/int(row_dict[ad_row[0]+ad_row[2]+ad_row[3]][i-2])
        else:
            row_dict[ad_row[0]+ad_row[2]+ad_row[3]]=ad_row

                    

init_row=['CROM','rsID','Ref','Alt','Type','Effect']
for i in prefs:
    init_row.extend(['AN_'+i,'AC_'+i,'AF_'+i])

with open ('ACE2_gnomad_all_extracted.tsv', 'w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter='\t') 
    my_writer.writerow(init_row)   
    for key in row_dict:
        my_writer.writerow(row_dict[key]) 
