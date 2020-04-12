# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 21:32:53 2018

@author: anton
"""
from collections import defaultdict
import csv
row_dict=defaultdict(list)

with open ('/home/anton/ACE2_gnomad_all_extracted.tsv', 'r',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t') 
    next(my_reader)  
    for row in my_reader:
        row_dict[row[0]+row[2]+row[3]]=row

eqtls_art=defaultdict(list)
prefs=['nfe','fin','eas','afr','sas','nfe_bgr','nfe_swe','nfe_nwe','nfe_est','nfe_seu','nfe_onf']
init_row=['CROM','rsID','Ref','Alt','Type','Effect']
for i in prefs:
    init_row.extend(['AN_'+i,'AC_'+i,'AF_'+i])
init_row.extend(['pval','tissue'])

with open ('/home/anton/eqtls_metadata.txt', 'r',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t') 
    next(my_reader)  
    for row in my_reader:
        eqtls_art[row[0]+row[2]+row[3]]=row[0:4]+[row[4].split(',')[0],row[4].split(',')[2]]


with open ('/home/anton/ACE2_eqtls.tsv', 'r',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')  
    for row in my_reader:
        key=row[0]+row[1]+row[2]
        if key  not in eqtls_art.keys():
            tissues=row[4].split('; ')
            pvals=[float(x) for x in row[5].split('; ')]
            #print(pvals[pvals.index(min(pvals))],tissues[pvals.index(min(pvals))])
            eqtls_art[key]=row[0:3]+[pvals[pvals.index(min(pvals))],tissues[pvals.index(min(pvals))].replace('_','-')]



with open ('/home/anton/ACE2_eqtls_gnomad.tsv', 'w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter='\t') 
    my_writer.writerow(init_row)
    for key in eqtls_art:
        if key in row_dict.keys():
            if len(eqtls_art[key])==5:
                my_writer.writerow(row_dict[key]+eqtls_art[key][3:5])
            else:
                my_writer.writerow(row_dict[key]+eqtls_art[key][4:6])
