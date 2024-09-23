"""
.. module:: extract_identity
    :platform: Any
    :synopsis: This module extracts information on the nucleotide in the 16 locations that determine identity
.. moduleauthor:: Paulo Nuin, August 2016
"""

import os
import subprocess
from pathlib import Path

from rich.console import Console

from .log_api import log_to_api

console = Console()


def mpileup(sample_id: str, datadir: str, identity: str, samtools: str) -> str:
    """
    Function that calls samtools to generate the pileup file that can be searched for the nucleotide state in
    the locations of interest.

    :param sample_id: ID of the patient/sample being analysed
    :param datadir: Location of the BAM files
    :param identity: Identity file name with defined location
    :param samtools: Samtools executable location

    :type sample_id: string
    :type datadir: string
    :type identity: string
    :type samtools: string

    :return: returns 'success' if the pileup file is successfully generated, 'exists' if the file already exists.

    :todo: return error
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"

    if Path(f"{bam_dir}/identity.mpileup").exists():
        console.log(f"Identity mpileup file exists {sample_id}")
        log_to_api(
            "Identity mpileup file exists",
            "INFO",
            "mpileup",
            sample_id,
            Path(datadir).name,
        )
        return "exists"

    console.log(f"Starting mpileup process for identity file {sample_id}")
    log_to_api(
        "Starting mpileup process for identity file",
        "INFO",
        "mpileup",
        sample_id,
        Path(datadir).name,
    )
    mpileup_string = f"{samtools} mpileup -l {identity} {bam_dir}{sample_id}.bam > {bam_dir}/identity.mpileup"
    console.log(mpileup_string)
    log_to_api(mpileup_string, "INFO", "mpileup", sample_id, Path(datadir).name)
    proc = subprocess.Popen(
        mpileup_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    console.log(f"mpileup generation completed {sample_id}")
    log_to_api(
        "mpileup generation completed", "INFO", "mpileup", sample_id, Path(datadir).name
    )
    return "success"


def get_nucleotides(pileup: str) -> tuple:
    """
    Function that returns nucleotide counts for the region

    :param pileup: reads under a defined location

    :type pileup: string

    :returns: returns the number of A, C, G and T
    """

    reads = "".join(e for e in pileup if e.isalnum()).upper()

    return reads.count("A"), reads.count("C"), reads.count("G"), reads.count("T")


def create_identity_table(sample_id: str, datadir: Path) -> None:
    """
    Function that creates the identity table in the sample data directory.

    :param sample_id: ID of the patient/sample being analysed
    :param datadir: Location where the table will be saved

    :type sample_id: str
    :type datadir: Path

    :return: None
    """

    mpileup_file = datadir / "BAM" / sample_id / "BAM" / "identity.mpileup"
    identity_file = datadir / "BAM" / sample_id / "identity.txt"

    with open(mpileup_file, "r") as f:
        mpileup = f.read().splitlines()

    if identity_file.exists():
        console.log(f"Identity file exists {sample_id}")
        log_to_api(
            "Identity file exists", "INFO", "mpileup", sample_id, datadir.name
        )
    else:
        console.log(f"Creating identity file {sample_id}")
        log_to_api(
            "Creating identity file", "INFO", "mpileup", sample_id, datadir.name
        )
        with open(identity_file, "w") as identity:
            for line in mpileup:
                temp = line.split()
                nucleotides = get_nucleotides(temp[4])
                identity.write(f"{temp[0]}\t{temp[1]}\t{str(nucleotides[0])}\t")
                identity.write(f"{str(nucleotides[1])}\t{str(nucleotides[2])}\t")
                identity.write(f"{str(nucleotides[3])}\n")
        console.log(f"Identity file creation complete {sample_id}")
        log_to_api(
            "Identity file creation complete",
            "INFO",
            "mpileup",
            sample_id,
            datadir.name
        )


if __name__ == "__main__":
    pass
