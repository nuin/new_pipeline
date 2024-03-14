"""
.. module:: variants_freebayes
    :platform: Any
    :synopsis: Module that generates variants by calling Freebayes
.. moduleauthor:: Paulo Nuin, January 2016

"""

import os
import subprocess
from pathlib import Path

from rich.console import Console

from .log_api import log_to_api

console = Console()



def freebayes_caller(datadir: str, sample_id: str, reference: str, bed_file: str, freebayes: str) -> str:
    """
    Function that calls Freebayes to generate a VCF file

    :param datadir: The directory where the data is located.
    :param sample_id: ID of the patient/sample being analysed
    :param reference: Reference file used in the original alignment
    :param bed_file: BED file with regions to be analysed
    :param freebayes: Location of Freebayes executable

    :type datadir: string
    :type sample_id: string
    :type reference: string
    :type bed_file: string
    :type freebayes: string

    :return: returns 'success' if the operation is successful, 'exists' if the VCF file already exists, 'error' if there is an error

    :rtype: string

    :todo: return error
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"
    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{vcf_dir}/{sample_id}_freebayes.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_freebayes.vcf file exists")
        log_to_api(
            "Freebayes VCF file exists",
            "info",
            "freebayes",
            sample_id,
            Path(datadir).name,
        )
        return "exists"

    console.log(f"Start variant calling with Freebayes {sample_id}")
    log_to_api(
        "Start variant calling with Freebayes",
        "info",
        "freebayes",
        sample_id,
        Path(datadir).name,
    )
    freebayes_string = (
        f"{freebayes} -f {reference} -v {vcf_dir}{sample_id}_freebayes.vcf -t {bed_file} -P 1 "
        f"{bam_dir}{sample_id}.bam"
    )
    console.log(f"Command {freebayes_string} {sample_id}")
    log_to_api(
        f"{freebayes_string}", "info", "freebayes", sample_id, Path(datadir).name
    )
    proc = subprocess.Popen(
        freebayes_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            console.log(output.decode("utf-8"))
    proc.wait()
    if os.path.getsize(f"{vcf_dir}{sample_id}_freebayes.vcf") == 0:
        console.log(f"Freebayes file size 0 {sample_id}", style="bold red")
        log_to_api(
            "Freebayes file size 0", "error", "freebayes", sample_id, Path(datadir).name
        )
        return "error"

    console.log(f"Freebayes variants determined {sample_id}")
    log_to_api(
        "Freebayes variants determined",
        "info",
        "freebayes",
        sample_id,
        Path(datadir).name,
    )
    return "success"


def edit_freebayes_vcf(sample_id: str, datadir: str) -> str:
    """
    Function that removes extra lines in Freebayes generated VCF to allow proper sorting.

    :param sample_id: ID of the patient/sample being analysed.
    :param datadir: The directory where the data is located.

    :type sample_id: string
    :type datadir: string

    :return: returns 'exists' if the VCF file already exists, 'success' if the operation is successful.

    :rtype: string
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"

    if Path(f"{vcf_dir}/{sample_id}_freebayes.final.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_freebayes.final.vcf file exists")
        log_to_api(
            "Freebayes VCF file exists",
            "info",
            "freebayes",
            sample_id,
            Path(datadir).name,
        )
        return "exists"

    freebayes_vcf = (
        open(f"{vcf_dir}/{sample_id}_freebayes.sorted.vcf").read().splitlines()
    )
    to_save = ""
    for line in freebayes_vcf:
        if not line.startswith("##contig=<ID"):
            to_save += line + "\n"

    console.log(f"Saving edited Freebayes VCF {sample_id}")
    log_to_api(
        "Saving edited Freebayes VCF",
        "info",
        "freebayes",
        sample_id,
        Path(datadir).name,
    )
    freebayes_final = open(f"{vcf_dir}/{sample_id}_freebayes.final.vcf", "w")
    freebayes_final.write(to_save)
    freebayes_final.close()
    console.log(f"Freebayes VCF edited {sample_id}")
    log_to_api(
        "Freebayes VCF edited", "info", "freebayes", sample_id, Path(datadir).name
    )

    return "success"
