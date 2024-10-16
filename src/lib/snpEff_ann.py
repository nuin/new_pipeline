"""
.. module:: snpEff_ann
    :platform: any
    :synopsis: This module calls snpEff to annotate variants from the merged VCF file
.. moduleauthor:: Paulo Nuin, June 2016

"""

import subprocess
from pathlib import Path

from rich.console import Console

from .log_api import log_to_api

console = Console()


def annotate_merged(sample_id: str, datadir: str, snpEff: str) -> str:
    """
    Function that calls snpEff to annotate the required VCF file.

    :param sample_id: ID of the patient/sample being analysed
    :param datadir: Location of the BAM files
    :param snpEff: snpEff executable location

    :type sample_id: string
    :type datadir: string
    :type snpEff: string

    :return: returns 'success' if the operation is successful, 'exists' if the VCF file already exists

    :rtype: string

    :todo: return error
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/{sample_id}"

    if Path(f"{vcf_dir}_merged.ann.vcf").exists():
        console.log(f"VCF file {vcf_dir}_merged.ann.vcf exists")
        log_to_api("VCF file exists", "INFO", "snpEff", sample_id, Path(datadir).name)
        return "exists"

    console.log(f"Starting snpEff annotation {sample_id}")
    log_to_api(
        "Starting snpEff annotation", "INFO", "snpEff", sample_id, Path(datadir).name
    )
    snpeff_string = f"{snpEff} hg19 {vcf_dir}_merged.vcf -onlyTr /apps/data/src/bundle/transcripts_only.txt"
    console.log(snpeff_string)
    log_to_api(snpeff_string, "INFO", "snpEff", sample_id, Path(datadir).name)
    proc = subprocess.Popen(
        snpeff_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    to_write = ""
    while True:
        output = proc.stdout.readline().strip()
        to_write += output.decode("utf-8") + "\n"
        if output == b"":
            break
        else:
            console.log(output.decode("utf-8")[:50])
    proc.wait()

    merged_ann = open(f"{vcf_dir}_merged.ann.vcf", "w")
    merged_ann.write(to_write)
    merged_ann.close()
    console.log("Annotated VCF file generated")
    log_to_api(
        "Annotated VCF file generated", "INFO", "snpEff", sample_id, Path(datadir).name
    )

    return "success"


if __name__ == "__main__":
    pass
