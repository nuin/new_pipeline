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
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_code(sample_id):

    return sample_id[-2:]


def octopus_caller(sample_id, directory, reference, bed_file, octopus):
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

    argument_bam = f"{directory}/BAM/{sample_id}/BAM/{sample_id}"
    argument_vcf = f"{directory}/BAM/{sample_id}/VCF/{sample_id}"

    if os.path.isfile(argument_vcf + "_octopus.vcf"):
        logger.info("Octopus VCF file exists " + sample_id)
    else:
        logger.info("Starting Octopus call " + sample_id)
        # platypus_string = '%s callVariants --refFile=%s --bamFiles=%s.recal_reads.bam --nCPU=8 --regions=%s -o %s_platypus_unsorted.vcf' % (
        # platypus, reference, argument, bed_file, argument3)
        octopus_string = "%s -R %s -I %s.recal_reads.bam -t %s -o %s_octopus.vcf" % (
            octopus,
            reference,
            argument_bam,
            bed_file,
            argument_vcf,
        )
        logger.info("Command " + octopus_string + " " + sample_id)
        proc = subprocess.Popen(
            octopus_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        while True:
            output = proc.stderr.readline().strip()
            if output == b"":
                break
            else:
                logger.info(output.decode("utf-8"))
        proc.wait()

        logger.info("Octopus VCF file created " + sample_id)
    change_vcf_version(argument_vcf + "_octopus.vcf")
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


# if __name__ == '__main__':

#     data_directory = '/Volumes/Jupiter/CancerPlusRuns/190207_NB551084_0058_AH2J7FAFXY_Cplus_2019_NGS_05_TEST'
#     sample_id = '19-017-020284C_KC_OS'
#     reference = '/opt/reference/hg19.fasta'
#     bed_file = '/opt/BED/Inherited_Cancer_panel_BED_91122_Target_adjusted_FINAL_GIPoly.bed'
#     octopus_ex = '/usr/local/bin/octopus'

#     octopus(sample_id, data_directory, reference, bed_file, octopus_ex)
