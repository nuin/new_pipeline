"""
.. module:: picard_actions
    :platform: any
    :synopsis: This module calls Picard to sort Freebayes VCF
.. moduleauthor:: Paulo Nuin, April 2016

"""

# pylint: disable-msg=too-many-arguments
# pylint: disable-msg=line-too-long

import subprocess
from pathlib import Path

from rich.console import Console

from .log_api import log_to_api

console = Console()


def picard_sort(datadir, sample_id, reference, picard):
    """
    Function that calls Picard for Freebayes VCF sorting

    :param sample_id: ID of the patient/sample being analysed
    :param datadir: Location of the BAM files
    :param reference: Reference file used in the original alignment
    :param picard: Locations of the picard jar file

    :type sample_id: string
    :type datadir: string
    :type reference: string
    :type picard: string
    """

    dictionary = reference.replace("fasta", "dict")

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"

    if Path(f"{vcf_dir}/{sample_id}_freebayes.final.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_freebayes.sorted.vcf file exists")
        log_to_api(
            "Freebayes sorted VCF file exists",
            "INFO",
            "picard_sort",
            sample_id,
            Path(datadir).name,
        )
        return "exists"

    console.log(f"Sorting Freebayes VCF result {sample_id}")
    log_to_api(
        "Sorting Freebayes VCF result",
        "INFO",
        "picard_sort",
        sample_id,
        Path(datadir).name,
    )
    picard_string = (
        f"{picard} SortVcf I={vcf_dir}/{sample_id}_freebayes.vcf "
        f"O={vcf_dir}/{sample_id}_freebayes.sorted.vcf "
        f"SEQUENCE_DICTIONARY={dictionary} QUIET=true"
    )
    console.log(f"Command {picard_string} {sample_id}")
    log_to_api(picard_string, "INFO", "picard_sort", sample_id, Path(datadir).name)
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            console.log(output)
    proc.wait()
    console.log(f"Freebayes sorted VCF file file created {sample_id}")
    log_to_api(
        "Freebayes sorted VCF file created",
        "INFO",
        "picard_sort",
        sample_id,
        Path(datadir).name,
    )

    return "success"


if __name__ == "__main__":
    datadir = "/Users/nuin/Projects/Data/Test_dataset"
    sample_id = "NA12877_1"
    reference = "/opt/reference/hg19.fasta"

    picard_sort(sample_id, datadir, reference, "picard")
