# Paulo Nuin July 2018
"""
To be added
"""
import sys
import pandas as pd


def print_full(x):
    """
    TBA
    """

    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')


def find_deletions(frame):
    """
    TBA
    """
    locations = pd.DataFrame()
    print('Deletions')
    for key, item  in frame.loc[frame['min'] <= 0.65].groupby('idmin'):
        locations = item
        print('Sample\t' + key)
        print('Gene', end='')
        [print('\t' + x) for x in set(item['gene'])]
        print('Windows\t' + str(len(item)))
        print('Unique Exons\t' + str(len(item.exon.unique())))
        [print(str(x) + '\t' + str(y) + '\tg.' + str(z) + '\tg.' + str(k) + '\t' + m) for x, y, z, k, m in zip(list(item['min']), list(item['var']), list(item['start']), list(item['end']), list(item['gene']))]
        print('\n')
        print(locations.groupby((locations.start.shift(-1) == locations.end).shift().ne(True).cumsum()).agg({'start': 'first', 'end': 'last', 'gene': 'last'}))

        print('\n')

def find_amplifications(frame):

    locations = pd.DataFrame()
    print('Amplifications')
    for key, item in frame.loc[frame['max'] >= 1.35].groupby('idmax'):
        locations = item
        print('Sample\t' + key)
        print('Gene', end='')
        [print('\t' + x) for x in set(item['gene'])]
        print('Windows\t' + str(len(item)))
        print('Unique Exons\t' + str(len(item.exon.unique())))
        [print(str(x) + '\t' + str(y) + '\tg.' + str(z) + '\tg.' + str(k) + '\t' + m) for x, y, z, k, m in zip(list(item['max']), list(item['var']), list(item['start']), list(item['end']), list(item['gene']))]
        print('\n')
        print(locations.groupby((locations.start.shift(-1) == locations.end).shift().ne(True).cumsum()).agg({'start': 'first', 'end': 'last', 'gene': 'last'}))


        print('\n')

def find_borderline(frame):

    print('Borderline')
    for key, item in frame.loc[(frame['min'] < 0.8) & (frame['min'] > 0.65)].groupby('idmin'):
        print('Sample\t' + key)
        print('Gene', end='')
        [print('\t' + x) for x in set(item['gene'])]
        print('\n')
        print('Windows\t' + str(len(item)))
        print('Unique Exons\t' + str(len(item.exon.unique())))
        print('Values:')
        [print(str(x) + '\t' + str(y) + '\tg.' + str(z) + '\tg.' + str(k) + '\t' + m) for x, y, z, k, m in zip(list(item['min']), list(item['var']), list(item['start']), list(item['end']), list(item['gene']))]

        print('\n')

    for key, item in frame.loc[(frame['min'] < 1.35) & (frame['min'] > 1.2)].groupby('idmax'):
        print('Sample\t' + key)
        print('Gene', end='')
        [print('\t' + x) for x in set(item['gene'])]
        print('\n')
        print('Windows\t' + str(len(item)))
        print('Unique Exons\t' + str(len(item.exon.unique())))
        print('Values:')
        [print(str(x) + '\t' + str(y) + '\tg.' + str(z) + '\tg.' + str(k) + '\t' + m) for x, y, z, k, m in zip(list(item['min']), list(item['var']), list(item['start']), list(item['end']), list(item['gene']))]

        print('\n')


def summarize(short_file, struct):
    """
    TBA
    """

    # my_data = open(short_file).read().splitlines()

    cnvs = pd.read_csv(short_file, sep='\t', index_col=0)

    cnvs['chr'] = [x[0] for x in list(cnvs.index.str.split(':'))]
    cnvs['gene'] = [x[1] for x in list(cnvs.index.str.split(':'))]
    cnvs['exon'] = [x[2].split('_')[0] for x in list(cnvs.index.str.split(':'))]
    cnvs['start'] = [x[2].split('_')[1] for x in list(cnvs.index.str.split(':'))]
    cnvs['end'] = [x[2].split('_')[2] for x in list(cnvs.index.str.split(':'))]


    if struct == 'deletion':
        find_deletions(cnvs)
    elif struct == 'amp':
        find_amplifications(cnvs)
    else:
        find_borderline(cnvs)


    # grouped = cnvs.groupby()


    # for i in cnvs.groupby('gene'):
    #     print(i)

    # cnvs['chr'] =

    # print(cnvs)

    # cnvs2 = cnvs.groupby((cnvs.start.shift(-1) == cnvs.end).shift().ne(True).cumsum())


if __name__ == '__main__':

    summarize(sys.argv[1], sys.argv[2])
