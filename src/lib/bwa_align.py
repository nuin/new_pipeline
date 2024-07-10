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
from pathlib import Path
from typing import List

from rich.console import Console

from .utils import move_bam

from .log_api import log_to_api

console = Console()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def run_bwa(
    sample_id: str,
    fastq_files: List[str],
    datadir: str,
    reference: str,
    bwa: str,
    samtools: str,
) -> int:
    """
    Performs the alignment/mapping for each pair/sample_id.

    Parameters:
    sample_id (str): The ID of the sample to be analysed.
    fastq_files (List[str]): The list of FASTQ files for the sample.
    datadir (str): The directory of the run.
    reference (str): The reference genome.
    bwa (str): The path to the BWA software.
    samtools (str): The path to the Samtools software.

    Returns:
    int: A flag indicating whether the BAM file exists (1) or was created (0).
    """

    console.log(f"Starting bwa processing for file {sample_id}")
    log_to_api(
        f"Starting bwa processing for file {sample_id}",
        "INFO",
        "bwa_align",
        sample_id,
        Path(datadir).name,
    )

    bam_index_check = 0
    if os.path.isfile(
        datadir + "/BAM/" + sample_id + "/" + sample_id + ".bam"
    ) or os.path.isfile(datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam"):
        console.log(f"{sample_id} BAM file exists")
        log_to_api(
            f"{sample_id} BAM file exists",
            "INFO",
            "bwa_align",
            sample_id,
            Path(datadir).name,
        )
        bam_index_check += 1
    elif os.path.isfile(
        datadir + "/BAM/" + sample_id + "/" + sample_id + ".bam.gz"
    ) or os.path.isfile(datadir + "/BAM/" + sample_id + "/BAM/" + sample_id + ".bam"):
        console.log(f"{sample_id} BAM file exists and compressed")
        log_to_api(
            f"{sample_id} BAM file exists and compressed",
            "INFO",
            "bwa_align",
            sample_id,
            Path(datadir).name,
        )
        bam_index_check += 1
    else:
        bam_header = f"@RG\\tID:{sample_id}\\tLB:{datadir}\\tPL:Illumina\\tSM:{sample_id}\\tPU:None"
        bwa_string = (
            f"{bwa} mem -t 16 -R '{bam_header}' {reference} {' '.join([str(file) for file in fastq_files])} "
            f"| {samtools} view -Sb - > {datadir}/BAM/{sample_id}/BAM/{sample_id}.bam"
        )
        console.log(f"{bwa_string}")
        log_to_api(
            f"{bwa_string}",
            "INFO",
            "bwa_align",
            sample_id,
            Path(datadir).name,
        )
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
        log_to_api(
            f"{samtools_string}",
            "INFO",
            "bwa_align",
            sample_id,
            Path(datadir).name,
        )
        console.log(f"Sorting BAM file {sample_id}")
        log_to_api(
            f"Sorting BAM file {sample_id}",
            "INFO",
            "bwa_align",
            sample_id,
            Path(datadir).name,
        )
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
        log_to_api(
            f"Indexing BAM file {sample_id}",
            "INFO",
            "bwa_align",
            sample_id,
            Path(datadir).name,
        )
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
        log_to_api(
            f"{sample_id} BAM file created",
            "INFO",
            "bwa_align",
            sample_id,
            Path(datadir).name,
        )

    return bam_index_check


if __name__ == "__main__":
    pass
