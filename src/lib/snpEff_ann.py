"""
.. module:: snpEff_ann
    :platform: any
    :synopsis: This module calls snpEff to annotate variants from the merged VCF file
.. moduleauthor:: Paulo Nuin, June 2016

"""

import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


def annotate_merged(sample_id, datadir, snpEff):
    """
    Function that calls snpEff to annotated the required
    VCF file

    :param sample_id: ID of the patient/sample being analysed
    :param datadir: Location of the BAM files
    :param snpEff: snpEff executable location

    :type sample_id: string
    :type datadir: string
    :type snpEff: string

    :return: returns success or exists

    :todo: return error
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/{sample_id}"

    if Path(f"{vcf_dir}_merged.ann.vcf").exists():
        console.log(f"VCF file {vcf_dir}_merged.ann.vcf exists")
        return "exists"

    console.log(f"Starting snpEff annotation {sample_id}")

    snpeff_string = (
        f"{snpEff} hg19 {vcf_dir}_merged.vcf -onlyTr /apps/data/src/bundle/transcripts_only.txt"
    )
    console.log(snpeff_string)

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

    return "success"


if __name__ == "__main__":

    datadir = "/Volumes/Juggernaut/200320_NB551084_0089_AH2YWNAFX2_Cplus_2020_NGS_07"
    # sample_id = '20GN-011G00002_PJ_DF'
    sample_id = "19-252-017170B_LL_FC"
    snpEff = "java -jar /usr/local/bin/snpEff.jar"
    # reference = '/opt/reference/hg19.fasta'
    # gatk = 'java -jar /usr/local/bin/GATK4.jar'

    # annotate_pseudo(sample_id, datadir, snpEff)
