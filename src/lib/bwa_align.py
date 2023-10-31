"""
.. module:: bwa_align
    :platform: Any
    :synopsis: Module that generates variants by calling Varscan
.. moduleauthor:: Paulo Nuin, July 2016
"""

import logging
import os
import subprocess
import time

from rich.console import Console

from lib.utils import move_bam

console = Console()

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

    console.log(f"Starting bwa processing for file {sample_id}")

    bam_index_check = 0
    if os.path.isfile(
        datadir + "/BAM/" + sample_id + "/" + sample_id + ".bam"
    ) or os.path.isfile(datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam"):
        console.log(f"{sample_id} BAM file exists")
        bam_index_check += 1
    elif os.path.isfile(
        datadir + "/BAM/" + sample_id + "/" + sample_id + ".bam.gz"
    ) or os.path.isfile(datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam"):
        console.log(f"{sample_id} BAM file exists and compressed")
        bam_index_check += 1
    else:
        bam_header = f"@RG\\tID:{sample_id}\\tLB:{datadir}\\tPL:Illumina\\tSM:{sample_id}\\tPU:None"
        bwa_string = (
            f"{bwa} mem -t 16 -R '{bam_header}' {reference} {' '.join(fastq_files)} "
            f"| {samtools} view -Sb - > {datadir}/BAM/{sample_id}/BAM/{sample_id}.bam"
        )
        console.log(f"{bwa_string}")
        proc = subprocess.Popen(
            bwa_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        while True:
            output = proc.stderr.readline().strip()
            if output == b"":
                break
            else:
                console.log(output.decode("utf-8"))
        proc.wait()
        time.sleep(5)
        samtools_string = (
            f"{samtools} sort -l 9 --threads 16 -o {datadir}/BAM/{sample_id}/BAM/{sample_id}.sorted.bam "
            f"{datadir}/BAM/{sample_id}/BAM/{sample_id}.bam"
        )
        console.log(f"{samtools_string}")
        console.log(f"Sorting BAM file {sample_id}")
        proc = subprocess.Popen(
            samtools_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
        )
        while True:
            output = proc.stderr.readline().strip()
            if output == b"":
                break
            else:
                console.log(output.decode("utf-8"))
        proc.wait()

        move_bam(datadir, sample_id, "sorted")

        console.log(f"Indexing BAM file {sample_id}")
        samtools_string_2 = (
            f"{samtools} index {datadir}/BAM/{sample_id}/BAM/{sample_id}.bam"
        )
        proc = subprocess.Popen(
            samtools_string_2,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )
        proc.wait()

        bam_index_check = 0
        console.log(f"{sample_id} BAM file created")

    return bam_index_check


if __name__ == "__main__":

    data_datadir = "/Users/nuin/Projects/Data/Test_dataset"
    sample_id = "NA12877_1"
    reference = "/opt/reference/hg19.fasta"
    run_bwa(
        sample_id,
        [
            "/Users/nuin/Projects/Data/Test_dataset/BaseCalls/NA12877_1_S1_L001_R1_001.fastq.gz",
            "/Users/nuin/Projects/Data/Test_dataset/BaseCalls/NA12877_1_S1_L001_R2_001.fastq.gz",
        ],
        data_datadir,
        reference,
        "bwa",
        "samtools",
    )
