"""
.. module:: variants_freebayes
    :platform: Any
    :synopsis: Module that generates variants by calling Freebayes
.. moduleauthor:: Paulo Nuin, January 2016

"""


import os
import subprocess
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_code(sample_id):

    return sample_id[-2:]


def freebayes_caller(sample_id, directory, reference, bed_file, freebayes):
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

    argument_bam = f"{directory}/BAM/{sample_id}/BAM/{sample_id}"
    argument_vcf = f"{directory}/BAM/{sample_id}/VCF/{sample_id}"

    if os.path.isfile(argument_vcf + "_freebayes.vcf"):
        logger.info("Freebayes VCF file exists " + sample_id)
        return "exists"

    logger.info("Starting Freebayes calling " + sample_id)
    freebayes_string = "%s -f %s -v %s_freebayes.vcf -t %s -P 1 %s.recal_reads.bam" % (
        freebayes,
        reference,
        argument_vcf,
        bed_file,
        argument_bam,
    )
    proc = subprocess.Popen(
        freebayes_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    logger.info("Command: " + freebayes_string + " " + sample_id)
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            logger.info(output.decode("utf-8") + " " + sample_id)
    proc.wait()
    if os.path.getsize(argument_vcf + "_freebayes.vcf") == 0:
        logger.error("Freebayes file size 0 " + sample_id)
        return "error"

    logger.info("INFO\t\tFrebayes call successful " + sample_id)
    return "success"


def edit_freebayes_vcf(sample_id, directory):
    """
    Function that removes extra lines in Freebayes generated VCF to allow proper sorting

    :param sample_id: ID of the patient/sample being analysed using GATK
    :param directory: Location of the BAM files

    :type sample_id: string
    :type directory: string
    """

    argument_vcf = f"{directory}/BAM/{sample_id}/VCF/{sample_id}"

    if os.path.isfile(argument_vcf + "_freebayes.final.vcf"):
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
