# Paulo Nuin Jan 2019

import sys
from collections import defaultdict
# import matplotlib.pyplot as plt
# from matplotlib_venn import venn3



result_table = open(sys.argv[1]).read().splitlines()

callers = defaultdict(dict)



for i in result_table:
    temp_dict = {}
    temp = i.split('\t')
    temp_dict['GATK'] = temp[4]
    temp_dict['FB'] = temp[5]
    temp_dict['PP'] = temp[5]
    callers[temp[0] + '_' + temp[1]] = temp_dict



all_callers = 0
g_f = 0
g_p = 0
f_p = 0
only_g = 0
only_f = 0
only_p = 0
no_one = 0
total = 0
for i in callers:
    total += 1
    if callers[i]['GATK'] == 'yes' and callers[i]['FB'] == 'yes' and callers[i]['PP'] == 'yes':
        all_callers += 1
    elif callers[i]['GATK'] == 'yes' and callers[i]['FB'] == 'yes' and callers[i]['PP'] == 'no':
        g_f += 1
    elif callers[i]['GATK'] == 'yes' and callers[i]['FB'] == 'no' and callers[i]['PP'] == 'yes':
        g_p += 1
    elif callers[i]['GATK'] == 'no' and callers[i]['FB'] == 'yes' and callers[i]['PP'] == 'yes':
        f_p += 1
    elif callers[i]['GATK'] == 'yes' and callers[i]['FB'] == 'no' and callers[i]['PP'] == 'no':
        only_g += 1
    elif callers[i]['GATK'] == 'no' and callers[i]['FB'] == 'yes' and callers[i]['PP'] == 'no':
        only_f += 1
    elif callers[i]['GATK'] == 'no' and callers[i]['FB'] == 'no' and callers[i]['PP'] == 'yes':
        only_p += 1
    else:
        no_one += 1


gatk = 0
total = 0
for i in callers:
    total += 1
    if callers[i]['GATK'] == 'yes':
        gatk += 1

fb = 0
total = 0
for i in callers:
    total += 1
    if callers[i]['FB'] == 'yes':
        fb += 1

pp = 0
total = 0
for i in callers:
    total += 1
    if callers[i]['PP'] == 'yes':
        pp += 1


print('GATK found: ' + str(gatk) + ' of 402')
print('Freebayes found: ' + str(fb) + ' of 402')
print('Platypus found: ' + str(pp) + ' of 402')
print('Found by all callers: ' + str(all_callers) + ' of 402')
print('Found by two callers: ' + str(g_f + g_p + f_p))
print('Found by GATK only: ' + str(only_g))
print('Found by Freebayes only: ' + str(only_f))
print('Found by platypus only: ' + str(only_p))

# venn3(subsets=(all_callers, g_f, g_p, f_p, no_one, 0, 0))
# venn3(subsets=(g_f, g_p, f_p, 6, 9, 4, all_callers))
# plt.show()

# print(total, all_callers, g_f, g_p, f_p, no_one)
        