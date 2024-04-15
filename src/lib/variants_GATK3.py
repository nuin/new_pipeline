"""
.. module:: variants_GATK3
    :platform: any
    :synopsis: Module that calls GATK to generate a VCF calls with variants and performs post-analysis of these variants
.. moduleauthor:: Paulo Nuin, January 2016

"""

import subprocess
from pathlib import Path

from rich.console import Console

from .log_api import log_to_api

console = Console()


def haplotype_caller(
    datadir: str, sample_id: str, reference: str, bed_file: str, gatk: str
) -> str:
    """
    Function that calls GATK's HaplotypeCaller parameter to generate a list of raw variants.

    :param datadir: The directory where the data is located.
    :param sample_id: ID of the patient/sample being analysed using GATK.
    :param reference: Reference file used in the original alignment.
    :param bed_file: BED file with regions to be analysed.
    :param gatk: GATK jar file location.

    :type datadir: string
    :type sample_id: string
    :type reference: string
    :type bed_file: string
    :type gatk: string

    :return: returns 'success' if the operation is successful, 'exists' if the VCF file already exists.

    :rtype: string
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"
    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{vcf_dir}/{sample_id}_GATK3.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_GATK.vcf file exists")
        log_to_api(
            "GATK VCF file exists", "info", "GATK", sample_id, Path(datadir).name
        )
        return "exists"

    console.log(f"Start variant calling with GATK3 {sample_id}")
    log_to_api(
        "Start variant calling with GATK3",
        "info",
        "GATK3",
        sample_id,
        Path(datadir).name,
    )
    GATK_string = (
        f"{gatk} -T HaplotypeCaller -R {reference} -I {bam_dir}/{sample_id}.bam -o {vcf_dir}{sample_id}_GATK3.vcf "
        f"-L {bed_file} -ip 2 -A StrandBiasBySample "
        f"-stand_call_conf 30 -pairHMM VECTOR_LOGLESS_CACHING -nct 16"
    )

    console.log(f"Command {GATK_string} {sample_id}")
    log_to_api(f"{GATK_string}", "info", "GATK3", sample_id, Path(datadir).name)
    proc = subprocess.Popen(
        GATK_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            console.log(output.decode("utf-8"))
    proc.wait()
    console.log("GATK3 variants determined " + sample_id)
    log_to_api(
        "GATK3 variants determined", "info", "GATK3", sample_id, Path(datadir).name
    )

    return "success"
