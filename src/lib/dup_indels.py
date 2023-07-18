"""
.. module:: dup_indels
    :platform: Any
    :synopsis: This module performs two call for Picard to remove duplicate reads and to add header information to the BAM files so they can be analysed in the remainder of the pipeline
.. moduleauthor:: Paulo Nuin, December 2015
"""

import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def remove_duplicates(sample_id, datadir, picard):
    """
    Function that runs the duplicate removal step in the analysis

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param datadir: Location of the BAM files
    :param picard: Picard jar file location

    :type sample_id: string
    :type datadir: string
    :type picard: string

    :return: returns success or exists

    """

    argument = datadir + '/BAM/' + sample_id + '/BAM/' + sample_id
    argument2 = datadir + '/BAM/' + sample_id + '/BAM/'

    if os.path.isfile(argument + '.dedup.bam') or os.path.isfile(argument2 + '.dedup.bam'):
        logger.info(argument + '.dedup.bam file exists ' + sample_id)
        return 'exists'

    logger.info('Picard - Starting duplicate removal')
    picard_string = '%s MarkDuplicates INPUT=%s.bam OUTPUT=%s.dedup.bam METRICS_FILE=%s.metrics.txt QUIET=true' % (picard, argument, argument, argument)
    logger.info('Command ' + picard_string + ' ' + sample_id)
    proc = subprocess.Popen(picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = proc.stderr.readline().strip()
        if output == b'':
            break
        else:
            logger.info(output.decode('utf-8'))
    proc.wait()
    if os.path.isfile(argument + '.dedup.bam'):
        logger.info('Picard - done duplicate marking ' + sample_id)
        return 'success'

    logger.error('dedup.bam, duplicate removal failed ' + sample_id)
    return 'error - process'


def samtools_duplicates(sample_id, datadir, samtools):
    """
    Function that runs the duplicate removal step in the analysis

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param datadir: Location of the BAM files
    :param samtools: Picard jar file location

    :type sample_id: string
    :type datadir: string
    :type samtools: string

    :return: returns success or exists

    """

    argument = datadir + '/BAM/' + sample_id + '/BAM/' + sample_id
    argument2 = datadir + '/BAM/' + sample_id + '/BAM/'

    try:
        os.remove(argument + '.pos.bam')
        logger.info('File removed ' + sample_id)
    except Exception as e:
        logger.info('File does not exist ' + str(e) + ' ' + sample_id)

    try:
        os.remove(argument + '.fix.bam')
        logger.info('File removed ' + sample_id)
    except Exception as e:
        logger.info('File does not exist ' + str(e) + ' ' + sample_id)

    try:
        os.remove(argument + '.sort.bam')
        logger.info('File removed ' + sample_id)
    except Exception as e:
        logger.info('File does not exist ' + str(e))

    if os.path.isfile(argument + '.dedup.bam') or os.path.isfile(argument2 + '.dedup.bam'):
        logger.info('Dedup.bam file exists ' + sample_id)
        return 'exists'

    if not os.path.isfile(argument + '.sort.bam'):
        logger.info('Samtools - Starting duplicate removal - sort 1')
        samtools_string = '%s sort -n -o %s.sort.bam %s.bam' % (samtools, argument, argument)
        logger.info('Command ' + samtools_string)
        proc = subprocess.Popen(samtools_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
    else:
        logger.info('Sort.bam file exists')

    if not os.path.isfile(argument + '.fix.bam'):
        logger.info('Samtools - Duplicate removal - fix')
        samtools_string2 = '%s fixmate -m %s.sort.bam %s.fix.bam' % (samtools, argument, argument)
        logger.info('Command ' + samtools_string2)
        proc = subprocess.Popen(samtools_string2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
    else:
        logger.info('Fix.bam file exists')

    if not os.path.isfile(argument + '.pos.bam'):
        logger.info('Samtools - Duplicate removal - sort 2')
        samtools_string3 = '%s sort -o %s.pos.bam %s.fix.bam' % (samtools, argument, argument)
        logger.info('Command ' + samtools_string3)
        proc = subprocess.Popen(samtools_string3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
    else:
        logger.info('Pos.bam file exists')

    if not os.path.isfile(argument + '.dedup.bam'):
        logger.info('Samtools - Duplicate removal - markdup')
        samtools_string4 = '%s markdup %s.pos.bam %s.dedup.bam' % (samtools, argument, argument)
        logger.info('Command ' + samtools_string4)
        proc = subprocess.Popen(samtools_string4, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
        return 'success'
    else:
        logger.info('Dedup.bam file exists')

    logger.error('dedup.bam, duplicate removal failed')
    return 'error - process'


def add_groups(sample_id, datadir, picard, samtools):
    """
    Function that add header groups to BAM files generated by BWA

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param datadir: Location of the BAM files
    :param samtools: Location of the Samtools executable

    :type sample_id: string
    :type datadir: string
    :type samtools: string

    :return: returns success or exists

    :todo: return error
    """

    if os.path.isdir(datadir + '/BAM/' + sample_id + '/BAM'):
        argument = datadir + '/BAM/' + sample_id + '/BAM/' + sample_id
        argument2 = datadir + '/BAM/' + sample_id + '/BAM/'
    else:
        argument = datadir + '/BAM/' + sample_id + '/' + sample_id
        argument2 = datadir + '/BAM/' + sample_id + '/'

    if os.path.isfile(argument + '.good.bam') or os.path.isfile(argument2 + '.good.bam'):
        logger.info('Good.bam file exists ' + sample_id)
        return 'exists'

    logger.info('Picard - adding groups to BAM file ' + sample_id)
    picard_string = '%s AddOrReplaceReadGroups INPUT=%s.dedup.bam OUTPUT=%s.good.bam RGSM=%s RGLB=Test RGPL=illumina RGPU=none QUIET=true' % (picard, argument, argument, sample_id)
    logger.info('Command ' + picard_string + ' ' + sample_id)
    samtools_string = '%s index %s.good.bam' % (samtools, argument)
    logger.info('Command ' + samtools_string + ' ' + sample_id)
    proc = subprocess.Popen(picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = proc.stderr.readline().strip()
        if output == b'':
            break
        else:
            logger.info(output.decode('utf-8'))
    proc.wait()

    if os.path.isfile(argument + '.good.bam'):
        logger.info('Picard - done adding groups ' + sample_id)
        logger.info('Starting samtools indexing of good.bam ' + sample_id)
        proc = subprocess.Popen(samtools_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
            output = proc.stderr.readline().strip()
            if output == b'':
                break
            else:
                logger.info(output.decode('utf-8'))
        proc.wait()
        return 'success'

    logger.info('adding groups failed ' + sample_id)
    return 'error - adding groups failed'


# if __name__ == '__main__':

#     data_datadir = '/Users/nuin/Projects/Data/Illumina_small/'
#     sample_id = 'NA12877_1'
#     remove_duplicates(sample_id, data_datadir, 'picard')
#     # add_groups(sample_id, data_datadir, 'picard', 'samtools')
