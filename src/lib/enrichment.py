"""
.. module:: extract_identity
    :platform: Any
    :synopsis: This module generates enrichment information
.. moduleauthor:: Paulo Nuin, August 2020
"""


import subprocess
from pathlib import Path

from dotenv import dotenv_values
from rich.console import Console

from .log_api import log_to_api

console = Console()


def get_enrichment(sample_id, datadir, panel):
    """

    :param sample_id:
    :param datadir:
    :return:
    """

    env = dotenv_values(f"{Path.cwd()}/.env")

    console.log(
        f"Saving enrichment file to {datadir}/BAM/{sample_id}/BAM/enrichment.enr"
    )
    log_to_api(
        f"Saving enrichment file to {datadir}/BAM/{sample_id}/BAM/enrichment.enr",
        "INFO",
        "extract_identity",
        sample_id,
        Path(datadir).name,
    )

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if panel == "Cplus":
        bed_file = "/apps/data/src/BED/new/C+_ALL_IDPE_01JUN2021.bed"
    else:
        bed_file = "/apps/data/src/BED/new/CardiacALL_29MAR2021.bed"

    enrichment_string = f"{env['ENRICHMENT']} {bam_dir}/{sample_id}.bam {bam_dir}/{sample_id}.bam.bai {bed_file} > {bam_dir}/enrichment.enr"

    if not Path(f"{bam_dir}/enrichment.enr").exists():
        console.log(
            f"Saving enrichment file to {datadir}/BAM/{sample_id}/BAM/enrichment.enr"
        )
        log_to_api(
            f"Saving enrichment file to {datadir}/BAM/{sample_id}/BAM/enrichment.enr",
            "INFO",
            "extract_identity",
            sample_id,
            Path(datadir).name,
        )
        console.log(enrichment_string)
        log_to_api(
            enrichment_string, "INFO", "extract_identity", sample_id, Path(datadir).name
        )
        proc = subprocess.Popen(
            enrichment_string,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        proc.wait()
        console.log(f"Enrichment file generated {sample_id}")
        log_to_api(
            "Enrichment file generated",
            "INFO",
            "extract_identity",
            sample_id,
            Path(datadir).name,
        )
    else:
        console.log(f"Enrichment file already exists {sample_id}")
        log_to_api(
            "Enrichment file already exists",
            "INFO",
            "extract_identity",
            sample_id,
            Path(datadir).name,
        )
