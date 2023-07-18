"""
.. module:: bwa_align
    :platform: Any
    :synopsis: Module that generates variants by calling Varscan
.. moduleauthor:: Paulo Nuin, July 2016
"""

import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def run_bwa(sample_id, fastq_files, datadir, reference, bwa, samtools):
    """
    function that performs the alignment/mapping, it is run iteratively for each pair/sample_id

    :param sample_id:
    :param fastq_files:
    :param datadir:
    :param reference:
    :param bwa:
    :param samtools:
    :return:
    """

    logger.info('Starting bwa processing for file ' + sample_id)

    bam_index_check = 0
    if os.path.isfile(datadir + "/BAM/" + sample_id + "/" + sample_id + ".bam") or os.path.isfile(datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam"):
        logger.info(sample_id + ' BAM file exists')
        bam_index_check += 1
    elif os.path.isfile(datadir + "/BAM/" + sample_id + "/" + sample_id + ".bam.gz") or os.path.isfile(datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam"):
        logger.info(sample_id + ' BAM file exists and compressed')
        bam_index_check += 1
    else:
        bwa_string = "%s mem -t 16 %s %s " % (bwa, reference, " ".join(fastq_files))
        bwa_string += "| %s view -Sb - | %s sort - -o %s" % (samtools, samtools, datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam")
        bwa_string += " && %s index %s" % (samtools, datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam")
        logger.info('Command ' + bwa_string)
        proc = subprocess.Popen(bwa_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
            output = proc.stderr.readline().strip()
            if output == b"":
                break
            else:
                logger.info(output.decode("utf-8"))
        # proc.wait()
        bam_index_check = 0
        logger.info('BWA alignment completed for sample ' + sample_id)

    return bam_index_check


if __name__ == '__main__':

    data_datadir = '/Users/nuin/Projects/Data/Test_dataset'
    sample_id = 'NA12877_1'
    reference = '/opt/reference/hg19.fasta'
    run_bwa(sample_id, ['/Users/nuin/Projects/Data/Test_dataset/BaseCalls/NA12877_1_S1_L001_R1_001.fastq.gz',
                         '/Users/nuin/Projects/Data/Test_dataset/BaseCalls/NA12877_1_S1_L001_R2_001.fastq.gz'], data_datadir, reference, 'bwa', 'samtools')