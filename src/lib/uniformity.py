"""
.. module:: uniformity
    :platform: any
    :synopsis: Generates uniformity values for samples
.. moduleauthor:: Paulo Nuin, Sept 2019

"""

import glob
import os
import sys
from pathlib import Path

import pandas as pd
from rich.console import Console

from .log_api import log_to_api

console = Console()


def get_coverage_values(datadir, panel):
    """ """

    if panel == "Cplus":
        panel_length = 328863
    elif panel == "Cardiac":
        panel_length = 725631

    console.log(f"Saving uniformity file to {datadir}/uniformity.txt")
    log_to_api(
        f"Saving uniformity file to {datadir}/uniformity.txt",
        "INFO",
        "uniformity",
        "NA",
        Path(datadir).name,
    )
    uniformity_file = open(f"{datadir}/uniformity.txt", "w")
    uniformity_file.write("Sample ID\t0.2*mean\t0.5*mean\t1.0*mean\tMean\n")
    for nucl_file in glob.glob(f"{datadir}/BAM/*/Metrics/*.nucl.out"):
        datadir = datadir.rstrip("/")
        sample_id = os.path.basename(nucl_file).replace(".nucl.out", "")
        values = pd.read_csv(nucl_file, sep="\t")
        cov_mean = values["coverage"].mean()
        cov_max = values["coverage"].max()
        values["cov_02"] = values["coverage"].between(cov_mean * 0.2, cov_max)
        values["cov_05"] = values["coverage"].between(cov_mean * 0.5, cov_max)
        values["cov_m"] = values["coverage"].between(cov_mean, cov_max)
        uniformity_file.write(sample_id + "\t")
        uniformity_file.write(
            str(values["cov_02"].value_counts().loc[True] / panel_length) + "\t"
        )
        uniformity_file.write(
            str(values["cov_05"].value_counts().loc[True] / panel_length) + "\t"
        )
        uniformity_file.write(
            str(values["cov_m"].value_counts().loc[True] / panel_length) + "\t"
        )
        uniformity_file.write(str(cov_mean) + "\n")
    uniformity_file.close()
    console.log(f"Uniformity file saved to {datadir}/uniformity.txt")
    log_to_api(
        f"Uniformity file saved to {datadir}/uniformity.txt",
        "INFO",
        "uniformity",
        "NA",
        Path(datadir).name,
    )


if __name__ == "__main__":

    get_coverage_values(sys.argv[1])
