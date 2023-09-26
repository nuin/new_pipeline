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

console = Console()


def get_code(sample_id):

    return sample_id[-2:]


def freebayes_caller(datadir, sample_id, reference, bed_file, freebayes):
    """
    Function that calls Freebayes to generate a VCF file

    :param sample_id: ID of the patient/sample being analysed
    :param directory: Location of the BAM files
    :param bed_file: BED file with regions to be analysed
    :param reference: Reference file used in the original alignment
    :param freebayes: Location of Freebayes executable

    :type sample_id: string
    :type directory: string
    :type bed_file: string
    :type reference: string
    :type freebayes: string

    :return: returns success or exists

    :todo: return error
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"
    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{vcf_dir}/{sample_id}_freebayes.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_freebayes.vcf file exists")
        return "exists"

    console.log(f"Start variant calling with Freebayes {sample_id}")
    freebayes_string = (
        f"{freebayes} -f {reference} -v {vcf_dir}{sample_id}_freebayes.vcf -t {bed_file} -P 1 "
        f"{bam_dir}{sample_id}.bam"
    )
    console.log(f"Command {freebayes_string} {sample_id}")
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
        return "error"

    console.log(f"Freebayes variants determined {sample_id}")
    return "success"


def edit_freebayes_vcf(sample_id, datadir):
    """
    Function that removes extra lines in Freebayes generated VCF to allow proper sorting

    :param sample_id: ID of the patient/sample being analysed using GATK
    :param directory: Location of the BAM files

    :type sample_id: string
    :type directory: string
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"

    if Path(f"{vcf_dir}/{sample_id}_freebayes.final.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_freebayes.final.vcf file exists")
        return "exists"

    freebayes_vcf = (
        open(f"{vcf_dir}/{sample_id}_freebayes.sorted.vcf").read().splitlines()
    )
    to_save = ""
    for line in freebayes_vcf:
        if not line.startswith("##contig=<ID"):
            to_save += line + "\n"

    console.log(f"Saving edited Freebayes VCF {sample_id}")
    freebayes_final = open(f"{vcf_dir}/{sample_id}_freebayes.final.vcf", "w")
    freebayes_final.write(to_save)
    freebayes_final.close()
    console.log(f"Freebayes VCF edited {sample_id}")

    return "success"
