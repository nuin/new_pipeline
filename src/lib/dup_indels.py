"""
.. module:: dup_indels
    :platform: Any
    :synopsis: This module performs two call for Picard to remove duplicate reads and to add header information to the BAM files so they can be analysed in the remainder of the pipeline
.. moduleauthor:: Paulo Nuin, December 2015
"""

import subprocess
from pathlib import Path

from rich.console import Console

from .log_api import log_to_api

console = Console()


def remove_duplicates(sample_id: str, datadir: str, picard: str) -> str:
    """
    Function that runs the duplicate removal step in the analysis using Picard.

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param datadir: Location of the BAM files
    :param picard: Picard jar file location

    :type sample_id: string
    :type datadir: string
    :type picard: string

    :return: returns 'success' if the process completes successfully, 'exists' if the dedup.bam file already exists, or 'error - process' if the duplicate removal fails.
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{bam_dir}/{sample_id}.dedup.bam").exists():
        console.log(f"{bam_dir}/{sample_id}.dedup.bam file exists")
        log_to_api(
            f"{bam_dir}/{sample_id}.dedup.bam file exists",
            "INFO",
            "dup_indels",
            sample_id,
            Path(datadir).name,
        )
        return "exists"
    elif Path(f"{bam_dir}/{sample_id}.dedup.bam.md5").exists():
        console.log(f"{bam_dir}/{sample_id}.dedup.bam.md5 file exists")
        log_to_api(
            f"{bam_dir}/{sample_id}.dedup.bam.md5 file exists",
            "INFO",
            "dup_indels",
            sample_id,
            Path(datadir).name,
        )
        return "exists"

    console.log("Picard - Starting duplicate removal")
    log_to_api(
        "Picard - Starting duplicate removal",
        "INFO",
        "dup_indels",
        sample_id,
        Path(datadir).name,
    )

    picard_string = (
        f"{picard} MarkDuplicates -I {bam_dir}/{sample_id}.bam -O {bam_dir}/{sample_id}.dedup.bam"
        f" -M {bam_dir}/{sample_id}.metrics.txt -MAX_RECORDS_IN_RAM 1000000 -AS true -QUIET true -CREATE_MD5_FILE true"
        f" -CREATE_INDEX true"
    )

    console.log(f"Command {picard_string} {sample_id}")
    log_to_api(
        f"Command {picard_string}", "INFO", "dup_indels", sample_id, Path(datadir).name
    )
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            console.log(output.decode("utf-8"))
    proc.wait()
    if Path(f"{bam_dir}/{sample_id}.dedup.bam").exists():
        console.log(f"Picard - done duplicate marking {sample_id}")
        log_to_api(
            "Picard - done duplicate marking",
            "INFO",
            "dup_indels",
            sample_id,
            Path(datadir).name,
        )
        return "success"

    console.log(f"dedup.bam, duplicate removal failed {sample_id}", style="bold red")
    log_to_api(
        "dedup.bam, duplicate removal failed",
        "ERROR",
        "dup_indels",
        sample_id,
        Path(datadir).name,
    )
    return "error - process"
