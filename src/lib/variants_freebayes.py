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


def edit_freebayes_vcf(sample_id, directory):
    """
    Function that removes extra lines in Freebayes generated VCF to allow proper sorting

    :param sample_id: ID of the patient/sample being analysed using GATK
    :param directory: Location of the BAM files

    :type sample_id: string
    :type directory: string
    """

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/"
    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{vcf_dir}/{sample_id}_freebayes.vcf").exists():
        console.log(f"{vcf_dir}/{sample_id}_GATK.vcf file exists")
        return "exists"

    logger.info("Editing sorted Freebayes VCF " + sample_id)
    freebayes_vcf = open(argument_vcf + "_freebayes.sorted.vcf").read().splitlines()
    to_save = ""
    for line in freebayes_vcf:
        if not line.startswith("##contig=<ID"):
            to_save += line + "\n"

    freebayes_final = open(argument_vcf + "_freebayes.final.vcf", "w")
    freebayes_final.write(to_save)
    freebayes_final.close()
    logger.info("Editing done " + sample_id)

    return "success"


# if __name__ == '__main__':

#     # data_directory = '/Users/nuin/Projects/Data/Test_dataset/'
#     # sample_id = 'NA12877_1'
#     # reference = '/opt/reference/hg19.fasta'
#     # bed_file = '/opt/BED/Inherited_Cancer_panel_FINAL.bed'

#     data_directory = '/Volumes/Jupiter/CancerPlusRuns/180531_NB551084_0040_AHW3TVAFXX_Cplus_2018_NGS_21/'
#     sample_id = '18-117-014682B_HL_DF'
#     reference = '/opt/reference/hg19.fasta'
#     bed_file = '/opt/BED/Inherited_Cancer_panel_FINAL.bed'

#     freebayes(sample_id, data_directory, reference, bed_file, 'freebayes')
#     edit_freebayes_vcf(sample_id, data_directory)
