# Paulo Nuin Sep 2023

import glob
from pathlib import Path
from shutil import move

from rich.console import Console

from .log_api import log_to_api

console = Console()


def move_bam(datadir: Path, sample: str, bam_file: str) -> bool:
    """
    Function that moves the BAM file and its index file from one location to another.

    :param datadir: The directory where the data is located.
    :param sample: The sample ID.
    :param bam_file: The name of the BAM file.

    :type datadir: Path
    :type sample: string
    :type bam_file: string

    :return: True if the operation is successful, False otherwise.

    :rtype: bool
    """

    bam_location = f"{datadir}/BAM/{sample}/BAM/{sample}"
    if Path(f"{bam_location}.{bam_file}.bam").exists():
        console.log(f"BAM file {bam_file}.bam exists")
        console.log(f"Moving {bam_location}.{bam_file}.bam to {bam_location}.bam")
        move(f"{bam_location}.{bam_file}.bam", f"{bam_location}.bam")
        try:
            move(f"{bam_location}.{bam_file}.bai", f"{bam_location}.bam.bai")
        except Exception as e:
            console.log(str(e))
            log_to_api(str(e), "ERROR", "move_bam", sample, datadir)
            console.log("Index file does not exist")
            log_to_api(
                "Index file does not exist", "ERROR", "move_bam", sample, datadir
            )
        console.log("Files moved")
        log_to_api("Files moved", "INFO", "move_bam", sample, datadir)
    else:
        console.log(
            f"BAM file {bam_file}.bam does not exist, but process might havbe been completed"
        )
        log_to_api(
            f"BAM file {bam_file}.bam does not exist, but process might havbe been completed",
            "WARNING",
            "move_bam",
            sample,
            datadir,
        )


def compile_identity(datadir: str) -> bool:
    """
    Function that reads samples' identity files and compiles them in a single file.

    :param datadir: The directory where the run is located.

    :type datadir: string

    :return: True if the operation is successful.

    :rtype: bool
    """

    all_identity = open(datadir + "/identity.txt", "w")
    for filename in glob.glob(datadir + "/BAM/*/*"):
        if filename.find("identity.txt") >= 0:
            all_identity.write(f"{filename.split('/')[-2]}\n")
            single_identity = open(filename).read()
            all_identity.write(single_identity)
    all_identity.close()

    return True
