from pathlib import Path
from shutil import move

from rich.console import Console

console = Console()


def is_bam_recalibrated(bam_path: Path) -> bool:
    """
    Check if the BAM file has already been recalibrated.
    This function can be implemented based on your specific criteria.
    """
    # For example, you could check for the existence of a recalibration log file
    recal_log = bam_path.with_name("recalibration.txt")
    return recal_log.exists()


def move_bam(datadir: Path, sample: str, bam_file: str) -> bool:
    bam_location = datadir / "BAM" / sample / "BAM"
    source_bam = bam_location / f"{sample}.{bam_file}.bam"
    target_bam = bam_location / f"{sample}.bam"
    source_bai = bam_location / f"{sample}.{bam_file}.bam.bai"
    target_bai = bam_location / f"{sample}.bai"

    if source_bam.exists():
        console.log(f"BAM file {bam_file}.bam exists")
        try:
            move(str(source_bam), str(target_bam))
            console.log(f"Moved {source_bam} to {target_bam}")

            if source_bai.exists():
                move(str(source_bai), str(target_bai))
                console.log(f"Moved {source_bai} to {target_bai}")
            else:
                console.log(f"Index file {source_bai} does not exist")

            return True
        except Exception as e:
            console.log(f"Error moving files: {str(e)}")
            return False
    else:
        console.log(f"BAM file {source_bam} does not exist")
        return False


def compile_identity(datadir: Path) -> bool:
    """
    Function that reads samples' identity files and compiles them in a single file.

    :param datadir: The directory where the run is located.

    :type datadir: Path

    :return: True if the operation is successful.

    :rtype: bool
    """

    all_identity_path = datadir / "identity.txt"
    with all_identity_path.open("w") as all_identity:
        for identity_file in datadir.glob("BAM/*/*/identity.txt"):
            sample_id = identity_file.parent.parent.name
            all_identity.write(f"{sample_id}\n")
            all_identity.write(identity_file.read_text())

    return True