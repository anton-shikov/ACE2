# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 21:32:53 2018

@author: anton
"""
from collections import defaultdict
import csv
import re


row_dict=defaultdict(list)
init_row=list()
with open ('/home/anton/ACE2_gnomad_all_extracted.tsv', 'r',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')   
    for row in my_reader:

        if row[0]=='POS':
            init_row=row
        else:
            row_dict[row[0]+row[2]+row[3]]=row
init_row.extend(['AN_ru','AC_rus','AF_rus'])

with open ('/home/anton/coronavirus_filtered_annotated.tsv', 'r',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')   
    next(my_reader)
    for row in my_reader:
        if row[2]+row[4]+row[5] in row_dict.keys():
            row_dict[row[2]+row[4]+row[5]].extend([row[23],row[22],row[24]])
        


with open ('ACE2_gnomad_rus.tsv', 'w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter='\t') 
    my_writer.writerow(init_row)   
    for key in row_dict:
        if len(row_dict[key]) > 39: 
            my_writer.writerow(row_dict[key]) 
