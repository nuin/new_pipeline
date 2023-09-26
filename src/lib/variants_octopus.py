"""
.. module:: variants_platypus
    :platform: Any
    :synopsis: Module that calls Octopus for variant calling
.. moduleauthor:: Paulo Nuin, February 2019

"""

# octopus -R /opt/reference/hg19.fasta -I
# /Volumes/Jupiter/CancerPlusRuns/190207_NB551084_0058_AH2J7FAFXY_Cplus_2019_NGS_05_TEST/BAM/19-015-021509B_SJ_OS/BAM/19-015-021509B_SJ_OS.recal_reads.bam
# -t /opt/BED/Inherited_Cancer_panel_BED_91122_Target_adjusted_FINAL_GIPoly.bed -o
# /Volumes/Jupiter/CancerPlusRuns/190207_NB551084_0058_AH2J7FAFXY_Cplus_2019_NGS_05_TEST/BAM/19-015-021509B_SJ_OS/VCF/19-015-021509B_SJ_OS.octopus.vcf

# pylint: disable-msg=too-many-arguments
# pylint: disable-msg=line-too-long

import os
import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


def get_code(sample_id):

    return sample_id[-2:]


def octopus_caller(datadir, sample_id, reference, bed_file, octopus):
    """
    Function that calls Platypus to generate a VCF file

    :param sample_id: ID of the patient/sample being analysed
    :param directory: Location of the BAM files
    :param bed_file: BED file with regions to be analysed
    :param reference: Reference file used in the original alignment
    :param octopus: Platypus Python script location

    :type sample_id: string
    :type directory: string
    :type bed_file: string
    :type reference: string
    :type platypus: string

    :return: returns success or exists

    :todo: return error
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"
    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{vcf_dir}/{sample_id}_octopus.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_octopus.vcf file exists")
        return "exists"

    console.log(f"Start variant calling with Octopus {sample_id}")
    octopus_string = (
        f"{octopus} -R {reference} -I {bam_dir}/{sample_id}.bam -t {bed_file} "
        f"-o {vcf_dir}/{sample_id}_octopus.vcf"
    )

    console.log(f"Command {octopus_string} {sample_id}")

    proc = subprocess.Popen(
        octopus_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            console.log(output.decode("utf-8"))
    proc.wait()

    console.log("Octopus VCF file created " + sample_id)
    change_vcf_version(f"{vcf_dir}{sample_id}_octopus.vcf")
    return "success"


def change_vcf_version(vcf_file):
    """
    Function that changes the VCF version from 4.3 to 4.2 so it can be used in GATK

    :param vcf_file: VCF file location

    :type vcf_file: string

    """

    contents = open(vcf_file).read().splitlines()

    contents[0] = contents[0].replace("4.3", "4.2")
    new_file = open(vcf_file, "w")
    for i in contents:
        new_file.write(i + "\n")
    new_file.close()
