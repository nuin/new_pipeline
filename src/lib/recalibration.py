"""
.. module:: recalibration
    :platform: Any
    :synopsis: Module that does base quality recalibration
.. moduleauthor:: Paulo Nuin, January 2016

"""

# pylint: disable-msg=too-many-arguments
# pylint: disable-msg=line-too-long

import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def base_recal1(sample_id, directory, bed_file, vcf_file, reference, gatk):
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

    argument = directory + '/BAM/' + sample_id + '/BAM/' + sample_id
    argument2 = directory + '/BAM/' + sample_id + '/BAM/'

    if os.path.isfile(argument2 + '/recal_data.table'):
        # log_to_kafka('Recalibration list exists ' + argument + ' ' + MY_HOSTNAME)
        logger.info('Recalibration list exists ' + sample_id)
        return 'exists'

    logger.info('Starting step one of base recalibration ' + sample_id)
    GATK_string = '%s BaseRecalibrator -R %s  -I %s.good.bam --known-sites %s -O %s/recal_data.table -L %s' % (gatk, reference, argument, vcf_file, argument2, bed_file)
    logger.info('Command ' + GATK_string +  ' ' + sample_id)
    proc = subprocess.Popen(GATK_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = proc.stderr.readline().strip()
        if output == b'':
            break
        else:
            logger.info(output.decode('utf-8') + ' ' + sample_id)
    proc.wait()
    logger.info('Recalibration step one completed successfully ' + sample_id)
    return 'success'


def base_recal2(sample_id, directory, vcf_file, reference, gatk):
    """
    Function that does the second step of the base recalibration
    process, generating the post-recalibration table

    :param sample_id: ID of the patient/sample being analysed using GATK
    :param directory: Location of the BAM files
    :param vcf_file: VCF file of known regions of variants
    :param reference: Reference file used in the original alignment
    :param gatk: GATK jar file location

    :type sample_id: string
    :type directory: string
    :type vcf_file: string
    :type reference: string
    :type gatk: string

    :todo: fix argument
    """

    argument = directory + '/BAM/' + sample_id + '/BAM/' + sample_id
    argument2 = directory + '/BAM/' + sample_id + '/BAM/'

    if os.path.isfile(argument2 + '/post_recal_data.table'):
        logger.info('Post recalibration list exists ' + sample_id)
        return 'exists'

    logger.info('Starting recalibration step two ' + sample_id)
    GATK_string = '%s ApplyBQSR -R %s -I %s.good.bam  %s/recal_data.table -O %s/post_recal_data.table' % (
        gatk, reference, argument, argument2, argument2)
    logger.info('Command ' + GATK_string + ' ' + sample_id)
    proc = subprocess.Popen(GATK_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = proc.stderr.readline().strip()
        if output == b'':
            break
        else:
            logger.info(output.decode('utf-8') + ' ' + sample_id)
    proc.wait()
    logger.info('Recalibration step two completed successfully ' + sample_id)
    return 'success'


def recalibrate(sample_id, directory, reference, gatk):
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
    argument = directory + '/BAM/' + sample_id + '/BAM/' + sample_id
    argument2 = directory + '/BAM/' + sample_id + '/BAM/'

    if os.path.isfile(argument + '.recal_reads.bam'):
        logger.info('Recalibrated BAM exists ' + sample_id)
        return 'exists'

    logger.info('Starting recalibration ' + sample_id)
    GATK_string = '%s ApplyBQSR -R %s -I %s.good.bam --bqsr-recal-file %s/recal_data.table -O %s.recal_reads.bam' % (
        gatk, reference, argument, argument2, argument)
    logger.info('Command ' + GATK_string + ' ' + sample_id)
    proc = subprocess.Popen(GATK_string, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = proc.stderr.readline().strip()
        if output == b'':
            break
        else:
            logger.info(output.decode('utf-8') + ' ' + sample_id)
    proc.wait()

    logger.info('Recalibration completed ' + sample_id)
    return 'success'


if __name__ == '__main__':

    data_directory = '/Users/nuin/Projects/Data/Test_dataset/'
    sample_id = 'NA12877_1'
    bed_file = '/opt/BED/Inherited_Cancer_panel_FINAL.bed'
    vcf_file = '/opt/bundle/dbsnp_138.hg19.vcf'
    reference = '/opt/reference/hg19.fasta'
    gatk = 'java -jar /usr/local/bin/GenomeAnalysisTK.jar'

    base_recal1(sample_id, data_directory, bed_file, vcf_file, reference, gatk)
    base_recal2(sample_id, data_directory, vcf_file, reference, gatk)
    recalibrate(sample_id, data_directory, reference, gatk)

# Later
# java -jar /usr/local/bin/GenomeAnalysisTK.jar -T AnalyzeCovariates -R  ~/New_projects/reference/hg19.fasta -before recal_data.table -after post_recal_data.table -csv report.csv

# Rscript ~/New_projects/bundle/BQSR.R report.csv recal_data.table test.pdf
