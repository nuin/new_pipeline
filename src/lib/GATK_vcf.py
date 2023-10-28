"""
.. module:: GATK_vcf
    :platform: any
    :synopsis: Module that calls GATK to compare VCF files and performs post-analysis of these variants
.. moduleauthor:: Paulo Nuin, April 2016

"""

import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


def vcf_comparison(datadir, sample_id, reference, gatk):
    """
    Function that merges the available VCFs in the sample VCF datadir

    :param sample_id: ID of the patient/sample being analysed using GATK
    :param datadir: Location of the BAM files
    :param bed_file: BED file with regions to be analysed
    :param reference: Reference file used in the original alignment
    :param gatk: GATK jar file location

    :type sample_id: string
    :type datadir: string
    :type bed_file: string
    :type reference: string
    :type gatk: string

    :return: returns success or exists
    """

    # --variant:varscan %svarscan_intersect.vcf
    # -genotypeMergeOptions PRIORITIZE -priority gatk,gatk3,freebayes,octopus

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF/{sample_id}"

    if Path(f"{vcf_dir}_merged.vcf").exists():
        console.log(f"VCF file {vcf_dir}_merged.vcf exists")
        return "exists"

    console.log(f"Starting merge GATK, Freebayes and Octopus VCFs {sample_id}")
    GATK_string = (
        f"{gatk} -T CombineVariants -R {reference} --variant:freebayes {vcf_dir}_freebayes.final.vcf"
        f" --variant:gatk {vcf_dir}_GATK.vcf --variant:gatk3 {vcf_dir}_GATK3.vcf "
        f"--variant:octopus {vcf_dir}_octopus.vcf -o {vcf_dir}_merged.vcf "
        f" --genotypemergeoption UNSORTED --mergeInfoWithMaxAC  --minimumN 2"
    )
    console.log(GATK_string)
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
    console.log("Merge GATK, Freebayes and Octopus VCFs: done")

    return "success"


# if __name__ == '__main__':

#     data_datadir = '/Users/nuin/Projects/Data/Test_dataset'
#     sample_id = 'NA12877_1'
#     vcf_comparison(sample_id, data_datadir, '/opt/reference/hg19.fasta',
#                    'java -jar /usr/local/bin/GenomeAnalysisTK.jar', '/Volumes/diagonalley/2017_Data/Cplus/NGS_test/')

#     # vcf_comparison('19-017-020284C_KC_OS', '/Volumes/Jupiter/CancerPlusRuns/190207_NB551084_0058_AH2J7FAFXY_Cplus_2019_NGS_05_TEST', '/opt/reference/hg19.fasta',
#     #                'java -jar /usr/local/bin/GenomeAnalysisTK.jar', '/Volumes/diagonalley/2017_Data/Cplus/NGS_test/')
