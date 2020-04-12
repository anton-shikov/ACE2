# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 19:48:05 2019

@author: anton
"""

import csv
import sys
from collections import defaultdict
from statistics import mean 

var_dict=defaultdict(list)
var_dict_ext=defaultdict(list)
init_row=[]

with open ('/home/anton/coronavirus_filtered.csv', 'r',newline='') as csvfile1:
    my_reader = csv.reader(csvfile1, delimiter=',')
    for row in my_reader:
        if 'Gene' in row:
            init_row=row[0:19]
            names_row=[row[i].replace('_GT','') for i in range(19,len(row),3)]
        else:
            mit_row=[]
            dep_row=[]
            for i in range(19,len(row),3):
                mit_row.append(row[i])
            for i in range(20,len(row),3):
                dep_row.append(row[i])
            dep_dat=[int(dep_row[i].split('/')[0])+int(dep_row[i].split('/')[1]) for i in range(len(dep_row)) if dep_row[i]!='./.']
            print(row[3],mean(dep_dat), dep_dat)
            key='|'.join([row[1],row[2],row[4],row[5]])
            #if mean(dep_dat)>=10:
            if True:
                #print(mean(dep_dat),row[3])
                var_dict[key]=row[0:19]+[mit_row.count('1/1'),mit_row.count('0/1'),mit_row.count('0/0'),mit_row.count('1/1')*2+mit_row.count('0/1'),mit_row.count('1/1')*2+mit_row.count('0/1')*2+mit_row.count('0/0')*2]+[mean(dep_dat)]
                var_dict_ext[key]=row[2:6]+[mit_row.count('1/1'),mit_row.count('0/1'),mit_row.count('0/0')]+[mean(dep_dat)]
                names11=[names_row[i] for i in range(len(mit_row)) if mit_row[i]=='1/1']
                names01=[names_row[i] for i in range(len(mit_row)) if mit_row[i]=='0/1']
            #depth=[names_row[i] for i in range(len(mit_row)) if mit_row[i]=='1/1']
                for i in range(len(names11)):
                    if names11[i][0]=='S':
                        names11[i]=names11[i][1:len(names11[i])]

                for i in range(len(names01)):
                    if names01[i][0]=='S':
                        names01[i]=names01[i][1:len(names01[i])]
                var_dict_ext[key].extend([', '.join(names11), ', '.join(names01),])


total_DP=list()
total_VAF=list()

with open ('/home/anton/coronavirus_filtered.Final.vcf', 'r',newline='') as csvfile1:
    my_reader = csv.reader(csvfile1, delimiter='\t')
    for row in my_reader:
        if '#' not in row[0]:
            key1='|'.join([row[0],row[1],row[3].split(':')[0],row[4].split(':')[0]])
            if key1 in var_dict.keys() and 'rs' in var_dict[key1][3]:
                MQ=row[7].split(';')[2].replace('MQ=','')
                if MQ=='NaN':
                    MQ=61.45
                anal_row=row[9:len(row)]
                DP=[int(anal_row[i].split(':')[2]) for i in range(len(anal_row)) if anal_row[i].split(':')[0]!='./.']
                DP_names=[var_dict[key1][3] for i in range(len(anal_row)) if anal_row[i].split(':')[0]!='./.']
                if len(DP)>1:
                   total_DP.append([DP,DP_names])
                VAF=[int(anal_row[i].split(':')[1].split(',')[1])/int(anal_row[i].split(':')[2]) for i in range(len(anal_row)) if anal_row[i].split(':')[0]=='0/1']
                VAF_names=[var_dict[key1][3] for i in range(len(anal_row)) if anal_row[i].split(':')[0]=='0/1']
                if len(VAF)>1:
                   total_VAF.append([VAF, VAF_names])

with open ('/home/anton/coronavirus_filtered_AF.frq', 'r',newline='') as csvfile1:
    my_reader = csv.reader(csvfile1, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        key1='|'.join([row[0],row[1],row[4].split(':')[0],row[5].split(':')[0]])
        if key1 in var_dict.keys():
            var_dict[key1].append(row[5].split(':')[1])
            var_dict_ext[key1].append(row[5].split(':')[1])


with open ('/home/anton/coronavirus_all_annotated.tsv', 'w',newline='') as csvfile1:
    my_reader = csv.writer(csvfile1, delimiter='\t')
    init_row.extend(['1/1 counts','0/1 counts', '0/0 counts','Allele Count','Allele Number','Mean_coverage','AF'])
    my_reader.writerow(init_row)
    for key in var_dict:
        my_reader.writerow(var_dict[key])


#with open ('/home/anton/coronavirus_filtered_samples.tsv', 'w',newline='') as csvfile1:
#    my_reader = csv.writer(csvfile1, delimiter='\t')
#    for key in var_dict_ext:
#        my_reader.writerow(var_dict_ext[key])


#with open ('/home/anton/VAF.tsv', 'w',newline='') as csvfile1:
#    my_reader = csv.writer(csvfile1, delimiter='\t')
#    for key in total_VAF:
#        for i in range(len(key[0])):
#            my_reader.writerow([key[0][i], key[1][i]])

#with open ('/home/anton/DP.tsv', 'w',newline='') as csvfile1:
#    my_reader = csv.writer(csvfile1, delimiter='\t')
#    for key in total_DP:
#        for i in range(len(key[0])):
#            my_reader.writerow([key[0][i], key[1][i]])


