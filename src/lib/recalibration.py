"""
.. module:: recalibration
    :platform: Any
    :synopsis: Module that does base quality recalibration
.. moduleauthor:: Paulo Nuin, January 2016

"""

# pylint: disable-msg=too-many-arguments
# pylint: disable-msg=line-too-long

import subprocess
from pathlib import Path

from rich.console import Console

from .log_api import log_to_api

console = Console()


def base_recal1(datadir, sample_id, bed_file, vcf_file, reference, gatk):
    """
    Function that does the first step of base recalibration, creating a
    recalibration data table

    :param sample_id: ID of the patient/sample being analysed
    :param directory: Location of the BAM files
    :param bed_file: BED file with regions to be analysed
    :param vcf_file: VCF file of known regions of variants
    :param reference: Reference file used in the original alignment
    :param gatk: GATK jar file location

    :type sample_id: string
    :type directory: string
    :type bed_file: string
    :type vcf_file: string
    :type reference: string
    :type gatk: string

    :return: returns success or exists

    :todo: return error
    :todo: fix argument
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{bam_dir}/recal_data.table").exists():
        console.log(f"{bam_dir}/recal_data.table file exists")
        log_to_api("recal_data.table file exists", "INFO", "base_recal1", sample_id, Path(datadir).name)
        return "exists"

    console.log("Starting step one of base recalibration for {sample_id}")
    log_to_api("Starting step one of base recalibration", "INFO", "base_recal1", sample_id, Path(datadir).name)

    GATK_string = (
        f"{gatk} BaseRecalibrator -R {reference}  "
        f"-I {bam_dir}/{sample_id}.bam --known-sites {vcf_file} "
        f"-O {bam_dir}/recal_data.table -L {bed_file}"
    )

    console.log(f"Command {GATK_string} {sample_id}")
    log_to_api(f"{GATK_string}", "INFO", "base_recal1", sample_id, Path(datadir).name)
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
    console.log(f"Recalibration step one completed successfully {sample_id}")
    log_to_api("Recalibration step one completed successfully", "INFO", "base_recal1", sample_id, Path(datadir).name)

    return "success"


def recalibrate(datadir, sample_id, reference, gatk):
    """
    Function that performs the third step of the base recalibration process,
    generating the final BAM file. *.recal_reads.bam

    :param sample_id: ID of the patient/sample being analysed using GATK
    :param directory: Location of the BAM files
    :param reference: Reference file used in the original alignment
    :param gatk: GATK jar file location

    :type sample_id: string
    :type directory: string
    :type reference: string
    :type gatk: string
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{bam_dir}/{sample_id}.recal_reads.bam").exists():
        console.log(f"{bam_dir}/{sample_id}.recal_reads.bam file exists")
        log_to_api("recal_reads.bam file exists", "INFO", "recalibrate", sample_id, Path(datadir).name)
        return "exists"
    elif Path(f"{bam_dir}/recalibration.txt").exists():
        console.log(f"{bam_dir}/recalibration.txt file exists")
        log_to_api("recalibration.txt file exists", "INFO", "recalibrate", sample_id, Path(datadir).name)
        return "exists"

    console.log(f"Starting recalibration {sample_id}")
    log_to_api("Starting recalibration", "INFO", "recalibrate", sample_id, Path(datadir).name)
    GATK_string = (
        f"{gatk} ApplyBQSR -R {reference} -I {bam_dir}/{sample_id}.bam "
        f"--bqsr-recal-file {bam_dir}/recal_data.table -O {bam_dir}/{sample_id}.recal_reads.bam"
    )
    console.log(f"Command {GATK_string} {sample_id}")
    log_to_api(f"{GATK_string}", "INFO", "recalibrate", sample_id, Path(datadir).name)
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

    console.log(f"Recalibration completed {sample_id}")
    log_to_api("Recalibration completed", "INFO", "recalibrate", sample_id, Path(datadir).name)

    recal_file = open(f"{bam_dir}/recalibration.txt", "w")
    recal_file.write(f"{GATK_string}\n")
    recal_file.close()

    return "success"
