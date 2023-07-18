"""
.. module:: variants_GATK3
    :platform: any
    :synopsis: Module that calls GATK to generate a VCF calls with variants and performs post-analysis of these variants
.. moduleauthor:: Paulo Nuin, January 2016

"""

import os
import subprocess
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_code(sample_id):

    return sample_id[-2:]


def haplotype_caller(sample_id, directory, reference, bed_file, gatk):
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

    argument_bam = f"{directory}/BAM/{sample_id}/BAM/{sample_id}"
    argument_vcf = f"{directory}/BAM/{sample_id}/VCF/{sample_id}"

    if os.path.isfile(argument_vcf + "_GATK3.vcf"):
        logger.info("Raw GATK3 variants file exists " + sample_id)
        return "exists"

    logger.info("Start variant calling with GATK3")
    GATK_string = (
        "%s -T HaplotypeCaller -R %s -I %s.recal_reads.bam  --genotyping_mode DISCOVERY -ip 10 "
        "-stand_call_conf 30 -o %s_GATK3.vcf -L %s -pairHMM VECTOR_LOGLESS_CACHING -nct 16 -A StrandBiasBySample"
        % (gatk, reference, argument_bam, argument_vcf, bed_file)
    )
    proc = subprocess.Popen(
        GATK_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    logger.info("Command:" + GATK_string + " " + sample_id)
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            logger.info(output.decode("utf-8"))
    proc.wait()
    logger.info("GATK3 variants determined " + sample_id)
    return "success"


# if __name__ == '__main__':

#     data_directory = '/Users/nuin/Projects/Data/Test_dataset/'
#     sample_id = 'NA12877_1'
#     reference = '/opt/reference/hg19.fasta'
#     bed_file = '/opt/BED/Inherited_Cancer_panel_FINAL.bed'
#     gatk = 'java -jar /usr/local/bin/GenomeAnalysisTK.jar'

#     haplotype_caller(sample_id, data_directory, reference, bed_file, gatk)
# do_CGP(sample_id, data_directory)
# filter_CGP(sample_id, data_directory)
# annotate(sample_id, data_directory)
