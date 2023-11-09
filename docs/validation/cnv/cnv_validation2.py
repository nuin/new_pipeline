# Paulo Nuin July 2018

import sys
import os
import pandas as pd


def print_full(x):
    """
    TBA
    """

    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')

def read_output(validation_file):
    """
    TBA
    """

    cnvs = pd.read_csv(validation_file, sep='\t', index_col=0, header=0)
    cnvs2 = pd.DataFrame(index=cnvs.index)

    cnvs2['var'] = cnvs.var(axis=1)
    cnvs2['max'] = cnvs.max(axis=1)
    cnvs2['min'] = cnvs.min(axis=1)

    print_full(cnvs.idxmax(axis=1))
    cnvs2['idmax'] = cnvs.idxmax(axis=1)
    cnvs2['idmin'] = cnvs.idxmin(axis=1)

    # cnvs['idmin'] = cnvs.idxmin(axis=1)
    output_file = os.path.basename(validation_file).replace('out', 'short')
    cnvs2.to_csv(output_file, sep='\t')

if __name__ == '__main__':

    
    print('process with sed -i \'\'  \'s/[[:blank:]]*$//\' $i before')
    read_output(sys.argv[1])

