# Paulo Nuin May 2018

import os
import sys
import glob

import pandas as pd
import vcf
from collections import defaultdict

check_v = defaultdict(dict)
variants_table = pd.read_csv(sys.argv[1], sep = '\t', header = 0)
variants_table = variants_table.dropna(axis = 0)
sample_table = open('sample_table').read().splitlines()
sample_dict = {}
for i in sample_table:
    temp = i.split('\t')
    sample_dict[temp[0]] = temp[1]

# variants_table['POS'] = variants_table['POS'].astype(int)
# for i,j in variants_table.itertuples(index = False):
#     print(a, b)

def get_vcf(sample_id, position, reference, alternate):

    bam_location = 'source/BAM'
    try:
        # platypus = vcf.Reader(open(bam_location + '/' + sample_dict[sample_id] + '/' + sample_dict[sample_id] + '_platypus.vcf'), 'r')
        # gatk = vcf.Reader(open(bam_location + '/' + sample_dict[sample_id] + '/' + sample_dict[sample_id] + '_GATK.vcf'), 'r')
        # freebayes = vcf.Reader(open(bam_location + '/' + sample_dict[sample_id] + '/' + sample_dict[sample_id] + '_freebayes.final.vcf'), 'r')
        merged = vcf.Reader(open(bam_location + '/' + sample_dict[sample_id] + '/' + sample_dict[sample_id] + '_merged.vcf'))
    except:
        print(bam_location + '/' + sample_dict[sample_id] + '/' + sample_dict[sample_id] + '_merged.vcf')
        print('file not found', sample_id)
        return 'file not found'

    # for item in merged:
    #     print(item.POS)
    #     print(item.INFO['set'])

    variants = {item.POS:item.INFO['set'] for item in merged}
    # variants_set = [item.INFO['set'] for item in merged]

    try:
        if variants[position] == 'Intersection':
            print(sample_id + '\t' + str(position) + '\t' + reference + '\t' + alternate + '\t' + 'yes\tyes\tyes')
        else:
            print(sample_id + '\t' + str(position) + '\t' + reference + '\t' + alternate + '\t' + variants[position])
    except:
        print(sample_id + '\t' + str(position) + '\t' + reference + '\t' + alternate + '\tnot found')


    # print(variants_pos, variants_set)


    #  if int(position) in variant_pos:
    #     print('found')

    # temp_dict = {}

    # variants_pos_g = [item.POS for item in gatk]
    # variants_pos_p = [item.POS for item in platypus]
    # variants_pos_f = [item.POS for item in freebayes]
    # if int(position) in variants_pos_g:
    #     temp_dict['GATK'] = 'found'
    #     check_v[sample_id + '_' + str(int(position))].update(temp_dict)
    #     print('found' + '\t' + sample_id + '\t' + str(int(position)))
    # else:
    #     temp_dict['GATK'] = 'NOT FOUND'
    #     check_v[sample_id + '_' + str(int(position))].update(temp_dict)
    #     print('NOT FOUND', + '\t' + sample_id + '\t' + str(int(position)))

    # if int(position) in variants_pos_p:
    #     temp_dict['Platypus'] = 'found'
    #     check_v[sample_id + '_' + str(int(position))].update(temp_dict)
    #     print('found' + '\t' + sample_id + '\t' + str(int(position)))
    # else:
    #     temp_dict['Platypus'] = 'NOT FOUND'
    #     check_v[sample_id + '_' + str(int(position))].update(temp_dict)
    #     print('NOT FOUND', + '\t' + sample_id + '\t' + str(int(position)))


    # if int(position) in variants_pos_f:
    #     temp_dict['Freebayes'] = 'found'
    #     check_v[sample_id + '_' + str(int(position))].update(temp_dict)
    #     print('found' + '\t' + sample_id + '\t' + str(int(position)))
    # else:
    #     temp_dict['Freebayes'] = 'NOT FOUND'
    #     check_v[sample_id + '_' + str(int(position))].update(temp_dict)
    #     print('NOT FOUND', + '\t' + sample_id + '\t' + str(int(position)))


    # for item in gatk:
    #     if int(position) == item.POS:
    #         return sample_id, int(position), reference, alternate, item.POS, item.REF, item.ALT
    #     else:
    #         return 'NOT FOUND', sample_id, int(position)
    # except Exception as e:
    #     return sample_id, int(position), 'NOT FOUND'

#     for bamdir in glob.glob('source/BAM/*'):
#         if sample_id in bamdir:
#             print(sample_id, bamdir)
#             sample = os.path.basename(bamdir)
#             # try:
#             gatk = vcf.Reader(open('source/BAM/' + sample + '/' + sample + '_GATK.vcf'), 'r')
#             # print(gatk)
#             for item in gatk:
#                 if int(position) == item.POS:
#                     print(item)
#                     return sample_id, int(position), reference, alternate, item.POS, item.REF, item.ALT
#             # except Exception as e:
#             #     return sample_id, int(position), 'NOT FOUND'

found = []
not_found = []
# print('Start')
for i, j, k, m in zip(variants_table['Sample'], variants_table['POS'], variants_table['REF'], variants_table['ALT']):
    result = get_vcf(i, j, k, m)

# for i in check_v:
#     print(i + '\t' + check_v['GATK'] + '\t' + check_v['Freebayes'] + '\t' + check_v['Platypus'])


#     try:
#         found.append(result)
#         print('found\t', end='')
#         for r in result:
#             print(str(r) + '\t', end='')
#     except Exception as e:
#         not_found.append((i, j, k, m, result))
#         print('NOT FOUND\t', end='')
#         for r in result:
#             print(str(r) + '\t', end='')
#     print()


# print('Found ' + str(len(found)))
# print('Not found ' + str(len(not_found)))
