"""
.. module:: variants_GATK3
    :platform: any
    :synopsis: Module that calls GATK to generate a VCF calls with variants and performs post-analysis of these variants
.. moduleauthor:: Paulo Nuin, January 2016

"""

import os
import subprocess
import logging
from pathlib import Path
from rich.console import Console

console = Console()


def haplotype_caller(datadir, sample_id, reference, bed_file, gatk):
    """
    Function that calls GATK's HaplotypeCaller parameter to generate a list
    of raw variants

    :param sample_id: ID of the patient/sample being analysed using GATK
    :param directory: Location of the BAM files
    :param bed_file: BED file with regions to be analysed
    :param reference: Reference file used in the original alignment
    :param gatk: GATK jar file location
    :param pcr: If True variant call is made on BWA original's BAM

    :type sample_id: string
    :type directory: string
    :type bed_file: string
    :type reference: string
    :type gatk: string
    :type pcr: bool

    :return: returns success or exists
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"
    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{vcf_dir}/{sample_id}_GATK3.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_GATK.vcf file exists")
        return "exists"

    console.log(f"Start variant calling with GATK3 {sample_id}")
    # GATK_string = (
    #     "%s -T HaplotypeCaller -R %s -I %s.recal_reads.bam  --genotyping_mode DISCOVERY -ip 10 "
    #     "-stand_call_conf 30 -o %s_GATK3.vcf -L %s -pairHMM VECTOR_LOGLESS_CACHING -nct 16 -A StrandBiasBySample"
    #     % (gatk, reference, argument_bam, argument_vcf, bed_file)
    # )

    GATK_string = (
        f"{gatk} -T HaplotypeCaller -R {reference} -I {bam_dir}/{sample_id}.bam -o {vcf_dir}{sample_id}_GATK3.vcf "
        f"-L {bed_file} -ip 2 -A StrandBiasBySample "
        f"-stand_call_conf 30 -pairHMM VECTOR_LOGLESS_CACHING -nct 16"
    )

    console.log(f"Command {GATK_string} {sample_id}")
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

    return "success"
