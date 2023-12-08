# Paulo Nuin June 2018

import sys
from collections import defaultdict
import pandas as pd


def print_rows(section, columns):

    for i, r in section.iterrows():
        print(r['Location'] + '\t', end='')
        for j in columns:
            print(str(r[j]), end='\t')
        print()


def print_full(x):
    """
    TBA
    """

    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')

def find_cnvs(cnv_file):
    """
    TBA
    """

    cnv_results = defaultdict(list)
    regions = {}
    cnvs = pd.read_csv(cnv_file, sep='\t')

    print('\t'.join(cnvs.columns.values))
    columns = list(cnvs)
    columns.remove('Location')

    for column in columns:
        under = cnvs.loc[cnvs[column] <= 0.65]
        over = cnvs.loc[cnvs[column] >= 1.35]
        borderline1 = cnvs.loc[(cnvs[column] > 1.2) & (cnvs[column] < 1.35)]
        borderline2 = cnvs.loc[(cnvs[column] > 0.65) & (cnvs[column] < 0.8)]

    print_rows(borderline1, columns)
    print_rows(borderline2, columns)
    print_rows(over, columns)
    print_rows(under, columns)


    # for i, r in over.iterrows():
    #     print(r['Location'] + '\t', end='')
    #     for j in columns:
    #         print(str(r[j]), end='\t')
    #     print()



    # if len(under[patient_id]) != 0:
    #     under_genes = find_genes(under[patient_id], genes)
    #     for gene in under_genes:
    #         gene_df = extract_gene(cnvs, gene)
    #         regions[gene] = get_nucleotides(
    #             gene_df[[patient_id, 'Location']], patient_id)
    #         spark = sparkline(gene_df.drop(['Gene', 'Location'], axis=1).dropna(axis=1), patient_id, directory, gene)
    #         if gene == 'TTN' or gene == 'DMD':
    #             big_sparkline(gene_df.drop(['Gene', 'Location'], axis=1).dropna(axis=1), patient_id, directory, gene)
    #         cnv_results[gene].append(spark)

    # if len(over[patient_id]) != 0:
    #     over_genes = find_genes(over[patient_id], genes)
    #     for gene in over_genes:
    #         gene_df = extract_gene(cnvs, gene)
    #         regions[gene] = get_nucleotides(
    #             gene_df[[patient_id, 'Location']], patient_id)
    #         spark = sparkline(gene_df.drop(['Gene', 'Location'], axis=1).dropna(axis=1), patient_id, directory, gene)
    #         if gene == 'TTN' or gene == 'DMD':
    #             big_sparkline(gene_df.drop(['Gene', 'Location'], axis=1).dropna(axis=1), patient_id, directory, gene)
    #         cnv_results[gene].append(spark)

    return cnv_results, regions


if __name__ == '__main__':


    cnv_file = sys.argv[1]
    find_cnvs(cnv_file)
