
import sys

html_vars = [ '_'.join(i.split('\t')[0:3:2]) for i in open(sys.argv[1]).read().splitlines()]

all_vars = ['_'.join(i.split('\t')[0:2]) for i in open(sys.argv[2]).read().splitlines()]



print(set(all_vars).difference(set(html_vars)))


# L86963_248525328', 'L84221_117122310', 'L92388_135438976', 'D69367_25671332', 'L76738_56204392