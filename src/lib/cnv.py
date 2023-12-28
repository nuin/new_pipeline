"""
.. module:: cnv.py
    :platform: any
    :synopsis: Module that calculates and generates normalized read information
.. moduleauthor:: Paulo Nuin, October 2017
"""

import glob
import os
import sys

import pandas as pd
import requests
import yaml
from rich.console import Console
from pathlib import Path
from .log_api import log_to_api

console = Console()

XLINKED = ["FANCB", "DMD", "EMD", "FHL1", "GLA", "LAMP2", "TAZ", "GPC3"]


def print_full(x):
    """

    :param x:
    :return:
    """

    pd.set_option("display.max_rows", len(x))
    print(x)
    pd.reset_option("display.max_rows")


def check_file_size(cnv_file, panel, expected_lines):
    """

    :param cnv_file:
    :param panel:
    :param expected_lines:
    :return:
    """

    file_lines = len(open(cnv_file).read().splitlines())
    if file_lines != expected_lines:
        return False

    return True


def compile_samples(datadir):
    """
    Function that checks if all .cnv files were generated and extracts the information in
    each file, returning a pandas DataFrame

    :param datadir: run location
    :type datadir: string

    :return pandas DataFrame with raw coverage data for each window
    """

    # if datadir.find('Cplus') >= 0:
    #     panel = 'cplus'
    #     expected_lines = 7068
    # else:
    #     expected_lines = 15677
    #     panel = 'cardiac'

    all_cnvs = pd.DataFrame()
    for sample in glob.glob(datadir + "/BAM/*"):
        console.log(f"Getting CNV information from {sample}")
        log_to_api(f"Getting CNV information from {sample}", "INFO", "CNV", "NA", Path(datadir).name)
        sample_id = os.path.basename(sample)
        # if not check_file_size(f"{sample}/{sample_id}.cnv", panel, expected_lines):
        cnv_sample = pd.read_csv(
            sample + "/" + sample_id + ".cnv", sep="\t", usecols=[sample_id]
        )
        location = pd.read_csv(
            sample + "/" + sample_id + ".cnv", sep="\t", usecols=["Location"]
        )

        all_cnvs[sample_id] = cnv_sample

    console.log("CNV information from all samples collected")
    log_to_api("CNV information from all samples collected", "INFO", "CNV", "NA", Path(datadir).name)

    all_cnvs["Location"] = location
    # saves file, just in case
    all_cnvs.to_csv(datadir + "/cnv_counts.txt", sep="\t", index=True, header=True)
    all_cnvs = all_cnvs.set_index(["Location"])

    return all_cnvs


# def scaler(cnvs):

#     scaler = preprocessing.StandardScaler()
#     scaled_df = scaler.fit_transform(cnvs)
#     scaled_df = pd.DataFrame(scaled_df)

#     return scaled_df


def cnv_calculation(datadir, cnvs, yaml_file):
    """
    Function that generates the CNV intra and inter-sample normalization

    :param datadir: run location
    :param cnvs: raw coverage data for each window
    :param yaml_file: run configuration file

    :type datadir: string
    :type cnvs: DataFrame
    :type yaml_file: string
    """

    # reads configuration
    console.log(f"Reading configuration file {yaml_file}")
    log_to_api(f"Reading configuration file {yaml_file}", "INFO", "CNV", "NA", Path(datadir).name)
    configuration = yaml.load(open(yaml_file).read(), Loader=yaml.FullLoader)
    # split the genders for X-linked calculation
    males, females = split_genders(configuration["Gender"])

    # remove some columns from the DataFrame
    console.log("Removing columns from the DataFrame")
    log_to_api("Removing columns from the DataFrame", "INFO", "CNV", "NA", Path(datadir).name)
    cnvs = cnvs.loc[:, ~cnvs.columns.str.contains("^Unnamed")]

    # intra-sample normalization
    console.log("Intra-sample normalization")
    log_to_api("Intra-sample normalization", "INFO", "CNV", "NA", Path(datadir).name)
    cnvs2 = cnvs.apply(lambda x: x / x.sum(), axis=0)
    # sca = scaler(cnvs2)
    # sca['Location'] = cnvs2.index
    # sca.set_index('Location')
    # sca.to_csv(datadir + '/zscoreS.txt', sep='\t')

    # save table, just in case
    cnvs2.to_csv(datadir + "/cnv_sum.txt", sep="\t")

    # inter-sample normalization
    console.log("Inter-sample normalization")
    log_to_api("Inter-sample normalization", "INFO", "CNV", "NA", Path(datadir).name)
    cnvs3 = cnvs2.apply(lambda x: x / x.mean(), axis=1)

    # calculating overall standard deviation
    console.log("Calculating standard deviation")
    log_to_api("Calculating standard deviation", "INFO", "CNV", "NA", Path(datadir).name)
    cnvs3["std"] = cnvs3.std(axis=1)
    # calculating female standard deviation
    console.log("Calculating female std deviation"))
    cnvs3["stdF"] = cnvs3[females].std(axis=1)
    # calculating male standard deviation
    console.log("Calculating male std deviation")
    cnvs3["stdM"] = cnvs3[males].std(axis=1)
    # saving final file
    console.log("Saving final file")
    log_to_api("Saving final file", "INFO", "CNV", "NA", Path(datadir).name)
    cnvs3.to_csv(f"{datadir}/cnv_mean.txt", sep="\t")

    # do CNV calculation for X-linked genes if there are more than 2 males in the run
    # todo: if there's only one male?
    if len(males) >= 2:
        cnvs_xlinked = cnvs.loc[cnvs.index.str.contains("|".join(XLINKED))]
        cnvs = cnvs[~cnvs.index.str.contains("FANCB")]
        cnvs_calculation_xlinked(datadir, cnvs_xlinked, cnvs2, yaml_file)

    return "done"


def split_genders(samples):
    """
    Function that returns a list of males and females for X-linked and other gender related checks

    :param samples: dictionary of samples and related gender

    :type samples: dictionary
    """

    # males, females = [], []
    # for i in samples:
    #     if list(i.values())[0] == 'Female':
    #         females.append(list(i.keys())[0])
    #     else:
    #         males.append(list(i.keys())[0])

    females = [k for k, v in samples.items() if v == "Female" or v == "F"]
    males = [k for k, v in samples.items() if v == "Male" or v == "M"]

    return males, females


def get_xlinked():
    """

    :return:
    """

    r = requests.get("http://10.106.108.24:8010/api/?chromosome=chrX")

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.json()


def cnvs_calculation_xlinked(datadir, cnvs_xlinked, cnv_sum, yaml_file):
    """
    Function that calculates CNVs for x-linked genes
    There are 4 DataFrames used in this function
        cnv and cnv2 and 3: normalization intra-gender
        cnv4: ormalization inter-gender

    :param datadir: run location
    :param cnvs_xlinked: columns for the X-linked genes
    :param cnv_sum: full data matrix, already summed
    :param yaml_file: run configuration file

    :type datadir: string
    :type cnvs_xlinked: DataFrame
    :type cnv_sum: DataFrame
    type yaml_file: string
    """

    match_string = []
    for gene in get_xlinked():
        match_string.append(gene["symbol"])

    to_match = "|".join(match_string)

    configuration = yaml.load(open(yaml_file).read(), Loader=yaml.FullLoader)

    males, females = split_genders(configuration["Gender"])

    cnv_males = cnvs_xlinked[males]
    cnv_females = cnvs_xlinked[females]

    cnv4 = cnv_sum
    cnv4 = cnv4.apply(lambda x: x / x.sum(), axis=0)

    cnv4["avgF"] = cnv4[females].mean(axis=1)
    cnv4["avgM"] = cnv4[males].mean(axis=1)

    cnv4[males] = cnv4[males].apply(lambda x: x / cnv4["avgF"])
    cnv4[females] = cnv4[females].apply(lambda x: x / cnv4["avgM"])

    cnvs_chrX = cnv4.loc[cnv4.index.str.contains(to_match)]
    cnvs_chrX.to_csv(datadir + "/cnv_mean_xlinked_cross.txt", sep="\t")

    cnv_males2 = cnv_males.apply(lambda x: x / x.sum(), axis=0)
    cnv_males3 = cnv_males2.apply(lambda x: x / x.mean(), axis=1)
    cnv_males3["stdM"] = cnv_males3.std(axis=1)

    cnv_females2 = cnv_females.apply(lambda x: x / x.sum(), axis=0)
    cnv_females3 = cnv_females2.apply(lambda x: x / x.mean(), axis=1)
    cnv_females3["stdF"] = cnv_females3.std(axis=1)

    cnvs_xlinked_comb = pd.concat([cnv_females3, cnv_males3], axis=1)

    cnvs_xlinked_comb.to_csv(datadir + "/cnv_mean_xlinked.txt", sep="\t")

    return "done"
