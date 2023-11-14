# Paulo Nuin Jan 2019

import sys


result_table = open(sys.argv[1]).read().splitlines()


temp = result_table[0].split('\t')
print('<table class="table table is-bordered is-striped is-narrow is-fullwidth">')
print('<thead>\n<tr>')
[print('<th>', i, '</th>') for i in temp]
print('</tr>\n</thead>\n<tbody>')

for i in range(1,len(result_table)):
    print('<tr>')
    temp = result_table[i].split('\t')
    [print('\t<td>', j, '</td>') for j in temp]
    print('</tr>')
print('</tbody>\n</table>')
