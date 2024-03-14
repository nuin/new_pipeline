"""
.. module:: cnv.py
    :platform: any
    :synopsis: Module that calculates and generates normalized read information
.. moduleauthor:: Paulo Nuin, October 2017
"""

import glob
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd
import requests
import yaml
from rich.console import Console

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


def check_file_size(cnv_file: str, panel: str, expected_lines: int) -> bool:
    """
    Checks if the number of lines in the cnv_file matches the expected_lines.

    Parameters:
    cnv_file (str): The path to the CNV file.
    panel (str): The panel being used. This parameter is not used in the function.
    expected_lines (int): The expected number of lines in the CNV file.

    Returns:
    bool: True if the number of lines in the CNV file matches the expected_lines, False otherwise.
    """

    # Get the number of lines in the CNV file
    file_lines = len(open(cnv_file).read().splitlines())

    # Check if the number of lines in the CNV file matches the expected_lines
    if file_lines != expected_lines:
        return False

    return True


import glob
from typing import Dict

import pandas as pd
from rich.console import Console

console = Console()


def compile_samples(datadir: str) -> pd.DataFrame:
    """
    Function that checks if all .cnv files were generated and extracts the information in
    each file, returning a pandas DataFrame

    Parameters:
    datadir (str): The directory of the run.

    Returns:
    pd.DataFrame: DataFrame with raw coverage data for each window.
    """

    # Initialize an empty DataFrame to store all CNV data
    all_cnvs = pd.DataFrame()

    # Iterate over all samples in the BAM directory
    for sample in glob.glob(datadir + "/BAM/*"):
        # Log the current sample being processed
        console.log(f"Getting CNV information from {sample}")

        # Extract the sample ID from the sample path
        sample_id = os.path.basename(sample)

        # Read the CNV data for the current sample
        cnv_sample = pd.read_csv(
            f"{sample}/{sample_id}.cnv", sep="\t", usecols=[sample_id]
        )

        # Read the location data for the current sample
        location = pd.read_csv(
            f"{sample}/{sample_id}.cnv", sep="\t", usecols=["Location"]
        )

        # Add the CNV data for the current sample to the DataFrame
        all_cnvs[sample_id] = cnv_sample

    # Log that all CNV information has been collected
    console.log("CNV information from all samples collected")

    # Add the location data to the DataFrame
    all_cnvs["Location"] = location

    # Set the index of the DataFrame to the location data
    all_cnvs = all_cnvs.set_index(["Location"])

    # Return the DataFrame with all CNV data
    return all_cnvs


# def scaler(cnvs):

#     scaler = preprocessing.StandardScaler()
#     scaled_df = scaler.fit_transform(cnvs)
#     scaled_df = pd.DataFrame(scaled_df)

#     return scaled_df


def cnv_calculation(datadir: str, cnvs: pd.DataFrame, yaml_file: str) -> str:
    """
    Function that generates the CNV intra and inter-sample normalization

    Parameters:
    datadir (str): The directory of the run.
    cnvs (pd.DataFrame): DataFrame with raw coverage data for each window.
    yaml_file (str): The run configuration file.

    Returns:
    str: A string indicating the completion of the function.
    """

    # reads configuration
    console.log(f"Reading configuration file {yaml_file}")
    log_to_api(
        f"Reading configuration file {yaml_file}",
        "INFO",
        "CNV",
        "NA",
        Path(datadir).name,
    )
    configuration = yaml.load(open(yaml_file).read(), Loader=yaml.FullLoader)
    # split the genders for X-linked calculation
    males, females = split_genders(configuration["Gender"])

    # remove some columns from the DataFrame
    console.log("Removing columns from the DataFrame")
    log_to_api(
        "Removing columns from the DataFrame", "INFO", "CNV", "NA", Path(datadir).name
    )
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
    log_to_api(
        "Calculating standard deviation", "INFO", "CNV", "NA", Path(datadir).name
    )
    cnvs3["std"] = cnvs3.std(axis=1)
    # calculating female standard deviation
    console.log("Calculating female std deviation")
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


def split_genders(samples: Dict[str, str]) -> Tuple[List[str], List[str]]:
    """
    Function that returns a list of males and females for X-linked and other gender related checks

    Parameters:
    samples (Dict[str, str]): A dictionary of samples and related gender.
                              The keys are sample names and the values are genders ("Male", "Female", "M", "F").

    Returns:
    Tuple[List[str], List[str]]: Two lists, the first one contains the names of male samples,
                                 the second one contains the names of female samples.
    """

    # Create a list of female samples
    females = [k for k, v in samples.items() if v == "Female" or v == "F"]

    # Create a list of male samples
    males = [k for k, v in samples.items() if v == "Male" or v == "M"]

    return males, females


def get_xlinked() -> Any:
    """
    Sends a GET request to a specific URL and returns the JSON response.

    Returns:
    Any: The JSON response from the GET request.
    """

    # Send a GET request to the specific URL
    r = requests.get("http://10.106.108.24:8010/api/?chromosome=chrX")

    # If the request was not successful, raise an exception and exit the program
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    # Return the JSON response from the GET request
    return r.json()


from typing import Any


def cnvs_calculation_xlinked(
    datadir: str, cnvs_xlinked: pd.DataFrame, cnv_sum: pd.DataFrame, yaml_file: str
) -> str:
    """
    Function that calculates CNVs for x-linked genes. There are 4 DataFrames used in this function
    cnv and cnv2 and 3: normalization intra-gender
    cnv4: normalization inter-gender

    Parameters:
    datadir (str): The directory of the run.
    cnvs_xlinked (pd.DataFrame): Columns for the X-linked genes.
    cnv_sum (pd.DataFrame): Full data matrix, already summed.
    yaml_file (str): The run configuration file.

    Returns:
    str: A string indicating the completion of the function.
    """

    # Create a list of gene symbols to match
    match_string = []
    for gene in get_xlinked():
        match_string.append(gene["symbol"])

    to_match = "|".join(match_string)

    # Load the configuration file
    configuration = yaml.load(open(yaml_file).read(), Loader=yaml.FullLoader)

    # Split the genders for X-linked calculation
    males, females = split_genders(configuration["Gender"])

    # Create DataFrames for male and female samples
    cnv_males = cnvs_xlinked[males]
    cnv_females = cnvs_xlinked[females]

    # Normalize the data intra-sample
    cnv4 = cnv_sum
    cnv4 = cnv4.apply(lambda x: x / x.sum(), axis=0)

    # Calculate the average for female and male samples
    cnv4["avgF"] = cnv4[females].mean(axis=1)
    cnv4["avgM"] = cnv4[males].mean(axis=1)

    # Normalize the data inter-sample
    cnv4[males] = cnv4[males].apply(lambda x: x / cnv4["avgF"])
    cnv4[females] = cnv4[females].apply(lambda x: x / cnv4["avgM"])

    # Filter the data for X-linked genes
    cnvs_chrX = cnv4.loc[cnv4.index.str.contains(to_match)]
    cnvs_chrX.to_csv(datadir + "/cnv_mean_xlinked_cross.txt", sep="\t")

    # Normalize the data for male samples
    cnv_males2 = cnv_males.apply(lambda x: x / x.sum(), axis=0)
    cnv_males3 = cnv_males2.apply(lambda x: x / x.mean(), axis=1)
    cnv_males3["stdM"] = cnv_males3.std(axis=1)

    # Normalize the data for female samples
    cnv_females2 = cnv_females.apply(lambda x: x / x.sum(), axis=0)
    cnv_females3 = cnv_females2.apply(lambda x: x / x.mean(), axis=1)
    cnv_females3["stdF"] = cnv_females3.std(axis=1)

    # Combine the normalized data for male and female samples
    cnvs_xlinked_comb = pd.concat([cnv_females3, cnv_males3], axis=1)

    # Save the combined data to a file
    cnvs_xlinked_comb.to_csv(datadir + "/cnv_mean_xlinked.txt", sep="\t")

    return "done"
