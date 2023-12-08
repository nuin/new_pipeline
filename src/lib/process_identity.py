"""
.. module:: process_identity
    :platform: any
    :synopsis: Module that processes identify files and return a full table with barcodes and other info
.. moduleauthor:: Paulo Nuin, March 2018

"""

# pylint: disable-msg=too-many-arguments
# pylint: disable-msg=line-too-long

import glob
import os

import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


# positions for all genes used
GENES = {
    29416481: "ALK_E29",
    215645464: "BARD1_E4",
    14187449: "XPC_E16",
    112164561: "APC_E16",
    55214348: "EGFR_E4",
    116435768: "MET_E20",
    116436022: "MET_E21",
    30999280: "WRN_E26",
    90967711: "NBN_E10",
    145742879: "RECQL4_E3",
    98209594: "PTCH1_E23",
    133208979: "POLE_E45",
    32936646: "BRCA2_E17",
    89858602: "FANCI_E37",
    89849480: "FANCA_E16",
    59763347: "BRIP1_E19",
}

# hgvs notation for all variants/position
HGVS = {
    29416481: "c.4472A>G",
    215645464: "c.1134G>C",
    14187449: "c.2815C>A",
    112164561: "c.1635G>A",
    55214348: "c.474C>T",
    116435768: "c.3912C>T",
    116436022: "c.4071G>A",
    30999280: "c.3222G>T",
    90967711: "c.1197T>C",
    145742879: "c.132A>G",
    98209594: "c.3944C>T",
    133208979: "c.6252A>G",
    32936646: "c.7806-14T>C",
    89858602: "c.3906T>C",
    89849480: "c.1501G>A",
    59763347: "c.2755T>C",
}

# codes used in the identity report
CODES = {
    "ALK_E29000T": 1,
    "ALK_E290C0T": 2,
    "ALK_E290C00": 3,
    "ALK_E29": 5,
    "BARD1_E40C00": 1,
    "BARD1_E40CG0": 2,
    "BARD1_E400G0": 3,
    "BARD1_E4": 5,
    "XPC_E1600G0": 1,
    "XPC_E1600GT": 2,
    "XPC_E16000T": 3,
    "XPC_E16": 5,
    "APC_E1600G0": 1,
    "APC_E16A0G0": 2,
    "APC_E16A000": 3,
    "APC_E16": 5,
    "EGFR_E40C00": 1,
    "EGFR_E40C0T": 2,
    "EGFR_E4000T": 3,
    "EGFR_E4": 5,
    "MET_E200C00": 1,
    "MET_E200C0T": 2,
    "MET_E20000T": 3,
    "MET_E20": 5,
    "MET_E2100G0": 1,
    "MET_E21A0G0": 2,
    "MET_E21A000": 3,
    "MET_E21": 5,
    "WRN_E2600G0": 1,
    "WRN_E2600GT": 2,
    "WRN_E26000T": 3,
    "WRN_E26": 5,
    "NBN_E10A000": 1,
    "NBN_E10A0G0": 2,
    "NBN_E1000G0": 3,
    "NBN_E10": 5,
    "RECQL4_E3000T": 1,
    "RECQL4_E30C0T": 2,
    "RECQL4_E30C00": 3,
    "RECQL4_E3": 5,
    "PTCH1_E2300G0": 1,
    "PTCH1_E23A0G0": 2,
    "PTCH1_E23A000": 3,
    "PTCH1_E23": 5,
    "POLE_E45000T": 1,
    "POLE_E450C0T": 2,
    "POLE_E450C00": 3,
    "POLE_E45": 5,
    "BRCA2_E17000T": 1,
    "BRCA2_E170C0T": 2,
    "BRCA2_E170C00": 3,
    "BRCA2_E17": 5,
    "FANCI_E37000T": 1,
    "FANCI_E370C0T": 2,
    "FANCI_E370C00": 3,
    "FANCI_E37": 5,
    "FANCA_E160C00": 1,
    "FANCA_E160C0T": 2,
    "FANCA_E16000T": 3,
    "FANCA_E16": 5,
    "BRIP1_E19A000": 1,
    "BRIP1_E19A0G0": 2,
    "BRIP1_E1900G0": 3,
    "BRIP1_E19": 5,

    "MET_E210CG0": "U",

}


def read_identity(sample_id, datadir):
    """
    Function that reads the full identity for each patient and returns a DataFrame with information

    :param sample_id: ID of the sample
    :param datadir: run location

    :type sample_id: string
    :type datadir: string

    :returns: DataFrame with sample data

    """

    sample_identity = pd.read_csv(
        f"{datadir}/BAM/{sample_id}/identity.txt",
        sep="\t",
        names=["Chromosome", "Position", "A", "C", "G", "T"],
    )

    return sample_identity


def generate_barcode(sample_id, datadir, sample_identity):
    """
    Function that reads the reads the codes and generates a numeric barcode for each sample

    :param sample_id: ID of the sample
    :param datadir: run_location
    :param sample_identity: list of nuclotides in each position of the barcode

    :type sample_id: string
    :type datadir: string


    """

    barcode = ""
    barcode_file = open(
        datadir + "/BAM/" + sample_id + "/" + sample_id + ".barcode.txt", "w"
    )
    for i in [0, 3, 1, 12, 15, 4, 14, 13, 5, 6, 8, 11, 10, 9, 7, 2]:
        barcode += str(CODES[sample_identity["Full code"][i]])

    barcode_file.write(barcode)
    console.log(f"Barcode for {sample_id} is {barcode}")


def process_identity(sample_id, datadir, sample_identity):
    """
    :todo: check function
    """

    sample_identity["ID"] = sample_identity["Position"].map(GENES)
    sample_identity["HGVS"] = sample_identity["Position"].map(HGVS)
    sample_identity["Total reads"] = sample_identity[["A", "C", "G", "T"]].sum(axis=1)
    sample_identity["pcA"] = sample_identity["A"] / sample_identity["Total reads"]
    sample_identity["pcC"] = sample_identity["C"] / sample_identity["Total reads"]
    sample_identity["pcG"] = sample_identity["G"] / sample_identity["Total reads"]
    sample_identity["pcT"] = sample_identity["T"] / sample_identity["Total reads"]
    sample_identity["codeA"] = np.where(sample_identity["pcA"] > 0.1, "A", "0")
    sample_identity["codeC"] = np.where(sample_identity["pcC"] > 0.1, "C", "0")
    sample_identity["codeG"] = np.where(sample_identity["pcG"] > 0.1, "G", "0")
    sample_identity["codeT"] = np.where(sample_identity["pcT"] > 0.1, "T", "0")
    sample_identity["Full code"] = (
        sample_identity["ID"]
        + sample_identity["codeA"]
        + sample_identity["codeC"]
        + sample_identity["codeG"]
        + sample_identity["codeT"]
    )

    sample_identity.to_csv(
        f"{datadir}/BAM/{sample_id}/{sample_id}.identity_full.txt",
        sep="\t",
        index=False,
    )
    console.log(f"Full identity file generated {sample_id}")

    return sample_identity


def barcoding(sample_id, datadir):
    """
    TBA
    """

    sample_identity = read_identity(sample_id, datadir)
    process_identity(sample_id, datadir, sample_identity)
    generate_barcode(sample_id, datadir, sample_identity)


def find_duplicates(datadir):
    """
    TBA
    """

    barcodes = open(datadir + "/barcodes.txt").read().splitlines()

    check = {}

    samples_dict = {}
    for line in barcodes:
        temp = line.split("\t")
        samples_dict[temp[0]] = temp[1]

    codes = list(samples_dict.values())
    for sample in samples_dict:
        if codes.count(samples_dict[sample]) > 1:
            check[sample] = True
        else:
            check[sample] = False

    return check


def compile_barcodes(datadir):
    """
    TBA
    """

    barcodes = open(datadir + "/barcodes.txt", "w")
    barcode_dict = {}
    for sample in glob.glob(datadir + "/BAM/*"):
        sample_id = os.path.basename(sample)
        barcode_file = (
            open(datadir + "/BAM/" + sample_id + "/" + sample_id + ".barcode.txt")
            .read()
            .splitlines()
        )
        barcode_dict[sample_id] = barcode_file[0]
        barcodes.write(sample_id + "\t" + barcode_file[0] + "\n")

    barcodes.close()

    return "success"


if __name__ == "__main__":

    # for sample_id in ['17-291-021140B_DM_OR',
    #                   '17-293-019731B_TL_DF',
    #                   '17-298-020323B_HM_DF',
    #                   '17-299-019925B_MJ_DF',
    #                   '18-268-019416B_CK_JA',
    #                   '18-268-020573B_RF_OS',
    #                   '18-268-020575B_WJ_48',
    #                   '18-269-020737B_WA_OS',
    #                   '18-269-020738B_BS_5Z',
    #                   '18-270-010956B_CJ_Z9',
    #                   '18-270-019017B_JM_DF',
    #                   '18-270-019024B_DJ_DF',
    #                   '18-270-020250B_HD_Z9',
    #                   '18-272-009307B_HJ_DF',
    #                   '18-283-010388B_MV_DF',
    #                   '18-283-010456B_LD_DF',
    #                   '18-283-010518B_HS_DF',
    #                   '18-285-006845B_M1_9P',
    #                   '18-285-006938B_M9_9P',
    #                   '18-285-007000B_M1_9P',
    #                   '18-285-008885B_M9_OR',
    #                   '18-285-008949B_M9_OR',
    #                   '18-285-009002B_M9_5Z',
    #                   '18-285-009062B_M9_OR',]:
    #   datadir = '/Volumes/Jupiter/CancerPlusRuns/181018_NB551084_0049_AHCH75AFXY_Cplus_2018_NGS_38/'
    #   barcoding(sample_id, datadir)
    datadir = "/Volumes/Jupiter/CancerPlusRuns/181018_NB551084_0049_AHCH75AFXY_Cplus_2018_NGS_38/"
    compile_barcodes(datadir)
