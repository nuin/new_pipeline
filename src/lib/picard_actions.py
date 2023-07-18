"""
.. module:: picard_actions
    :platform: any
    :synopsis: This module calls Picard to sort Freebayes VCF
.. moduleauthor:: Paulo Nuin, April 2016

"""

# pylint: disable-msg=too-many-arguments
# pylint: disable-msg=line-too-long

import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_code(sample_id):

    return sample_id[-2:]


def picard_sort(sample_id, datadir, reference, picard):
    """
    Function that calls Picard for Freebayes VCF sorting

    :param sample_id: ID of the patient/sample being analysed
    :param datadir: Location of the BAM files
    :param reference: Reference file used in the original alignment
    :param picard: Locations of the picard jar file

    :type sample_id: string
    :type datadir: string
    :type reference: string
    :type picard: string
    """

    dictionary = reference.replace('fasta', 'dict')

    code = get_code(sample_id)
    argument_vcf = f"{datadir }/BAM/{sample_id}/VCF/{sample_id}"

    if os.path.isfile(argument_vcf + '_freebayes.sorted.vcf'):
        logger.info('Freebayes sorted VCF file exists ' + sample_id)
        return 'exists'

    logger.info('Sorting Freebayes VCF result')
    picard_string = '%s SortVcf I=%s_freebayes.vcf O=%s_freebayes.sorted.vcf SEQUENCE_DICTIONARY=%s QUIET=true' % (picard, argument_vcf, argument_vcf, dictionary)
    logger.info(picard_string + ' ' + sample_id)
    proc = subprocess.Popen(picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = proc.stderr.readline().strip()
        if output == b'':
            break
        else:
            logger.info(output.decode('utf-8'))
    proc.wait()
    logger.info('Freebayes sorted VCF file file created ' + sample_id)
    return 'success'


if __name__ == '__main__':

    datadir = '/Users/nuin/Projects/Data/Test_dataset'
    sample_id = 'NA12877_1'
    reference = '/opt/reference/hg19.fasta'

    picard_sort(sample_id, datadir, reference, 'picard')
