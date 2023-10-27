"""
.. module:: extract_identity
    :platform: Any
    :synopsis: This module generates enrichment information
.. moduleauthor:: Paulo Nuin, August 2020
"""


import glob
import os
import subprocess
import sys
from pathlib import Path

from rich.console import Console

from dotenv import dotenv_values

console = Console()

def get_enrichment(sample_id, datadir, panel):
    """

    :param sample_id:
    :param datadir:
    :return:
    """

    env = dotenv_values(f"{Path.cwd()}/.env")

    console.log(f"Saving enrichment file to {datadir}/BAM/{sample_id}/BAM/enrichment.enr")

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if panel == "Cplus":
        bed_file = "/apps/data/src/BED/new/C+_ALL_IDPE_01JUN2021.bed"
    else:
        bed_file = "/apps/data/src/BED/new/CardiacALL_29MAR2021.bed"

    enrichment_string = f"{env['ENRICHMENT']} {bam_dir}/{sample_id}.bam {bam_dir}/{sample_id}.bam.bai {bed_file} > {bam_dir}/enrichment.enr"


    if not Path(f"{bam_dir}/enrichment.enr").exists():
        console.log(f"Saving enrichment file to {datadir}/BAM/{sample_id}/BAM/enrichment.enr")
        console.log(enrichment_string)
        proc = subprocess.Popen(enrichment_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
        console.log(f"Enrichment file generated {sample_id}")
    else:
        console.log(f"Enrichment file already exists {sample_id}")

