#!/usr/bin/env python3.6

import gzip
import sys
import re
import numpy as np
import math

populations = ['nfe', 'fin', 'eas', 'afr', 'sas', 'nfe_bgr', 'nfe_swe', 'nfe_est', 'nfe_seu', 'nfe_onf', 'nfe_nwe']

exomes, genomes, selected, perm_count, outfile = sys.argv[1:]
perm_count = int(perm_count)

variants = {}
var_info = {}

snp_counter = 0

with gzip.open(exomes, 'rt') as ex_hand:
    for line in ex_hand:
        if line.startswith('#'):
            continue
        snp_counter += 1
#        if snp_counter % 10000 == 0:
#            print(f'Processed {snp_counter} SNPs...')
#        if snp_counter > 300000:
#            finita = True
#            break
        content = line.strip().split('\t')
        var_id = f'{content[0]}:{content[1]}:{content[3]}:{content[4]}'
        if 'AC_nfe=' not in line or 'AN_nfe=' not in line:
            continue
        ac = int(re.findall('AC_nfe=(\d+)', content[7])[0])
        an = int(re.findall('AN_nfe=(\d+)', content[7])[0])
        if ac == 0:
            continue
        variants[var_id] = (ac, an)
        var_info[var_id] = dict()
        for pop in populations:
            ac = int(re.findall(f'AC_{pop}=(\d+)', content[7])[0]) if f'AC_{pop}' in line else 0
            an = int(re.findall(f'AN_{pop}=(\d+)', content[7])[0]) if f'AN_{pop}' in line else 0
            var_info[var_id][pop] = (ac, an)

#print(snp_counter)

with gzip.open(genomes, 'rt') as gen_hand:
    for line in gen_hand:
#        if finita:
#            break
        if line.startswith('#'):
            continue
        snp_counter += 1
#        if snp_counter % 10000 == 0:
#            print(f'Processed {snp_counter} SNPs...')
        content = line.strip().split('\t')
        var_id = f'{content[0]}:{content[1]}:{content[3]}:{content[4]}'
        if 'AC_nfe=' not in line or 'AN_nfe=' not in line:
            continue
        ac = int(re.findall('AC_nfe=(\d+)', content[7])[0])
        if ac == 0:
            continue
        an = int(re.findall('AN_nfe=(\d+)', content[7])[0])
        if var_id in variants:
            var_data = variants[var_id]
            var_data = (var_data[0] + ac, var_data[1] + an)
            variants[var_id] = var_data
        else:
            variants[var_id] = (ac, an)
        if var_id not in var_info:
            var_info[var_id] = dict()
        for pop in populations:
            ac = int(re.findall(f'AC_{pop}=(\d+)', content[7])[0]) if f'AC_{pop}' in line else 0
            an = int(re.findall(f'AN_{pop}=(\d+)', content[7])[0]) if f'AN_{pop}' in line else 0
            if pop in var_info[var_id]:
                var_info[var_id][pop] = (var_info[var_id][pop][0] + ac, var_info[var_id][pop][1] + an)
            else:
                var_info[var_id][pop] = (ac, an)

ace_vars = []
with open(selected, 'r') as sel_handle:
    for line in sel_handle:
        if line.startswith('POS'):
            continue
        content = line.strip().split('\t')
        ac = int(content[7])
        an = int(content[6])
        if ac == 0:
            continue
#        if af_bin == 0:
#            continue
        ace_vars.append((ac, an))

#print(freq_counts)
#print(freq_counts)

gnomad_binned = {}
for var_id in variants:
    var_data = variants[var_id]
    if var_data[1] == 0 or var_data[0] == 0:
        continue
#    if var_af_bin == 0:
#        continue
    gnomad_binned[var_data[0]] = gnomad_binned.get(var_data[0], []) + [var_id]

matchup = {}
for i in range(len(ace_vars)):
    matchup[i] = []
    ref_an = ace_vars[i][1]
    similar_ac = gnomad_binned[ace_vars[i][0]]
    for var_id in similar_ac:
        var_data = variants[var_id]
        if var_data[1] in range(int(ref_an - (ref_an/20)), int(ref_an + ref_an/20) + 1):
            matchup[i].append(var_id)

header_string = ''
for pop in populations:
    if 'nfe' in pop:
        header_string += f'{pop}\t'

random_file = open(outfile, 'w')

print(header_string.strip(), file=random_file)

for i in range(perm_count):
    this_iter = {pop: 0 for pop in populations if 'nfe' in pop}
    out_string = ''
    snp_counter = 0
    for snp in range(len(ace_vars)):
#        print(fbin)
        var_id = np.random.choice(matchup[i], size=1)[0]
        snp_counter += 1
        for pop in this_iter:
            pop_data = var_info[var_id][pop]
            if pop_data[1] > 0:
                this_iter[pop] += pop_data[0]/pop_data[1]
    for pop in this_iter:
        out_string += f'{this_iter[pop]}\t'
    if out_string:
        print(out_string.strip(), file=random_file)
    print(snp_counter)

print(var_info[var_id]['nfe'][0] == variants[var_id][0])


