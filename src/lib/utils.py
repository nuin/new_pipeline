# Paulo Nuin Sep 2023

from pathlib import Path
from rich.console import Console
from shutil import move

console = Console()


def move_bam(datadir: Path, sample: str, bam_file: str) -> bool:
    """

    :param datadir:
    :param bam_file:
    :return:
    """

    bam_location = f"{datadir}/BAM/{sample}/BAM/{sample}"
    if Path(f"{bam_location}.{bam_file}.bam").exists():
        console.log(f"BAM file {bam_file}.bam exists")
        console.log(f"Moving {bam_location}.{bam_file}.bam to {bam_location}.bam")
        move(f"{bam_location}.{bam_file}.bam", f"{bam_location}.bam")
        move(f"{bam_location}.{bam_file}.bai", f"{bam_location}.bam.bai")
        console.log("Files moved")
    else:
        console.log(f"BAM file {bam_file}.bam does not exist")
