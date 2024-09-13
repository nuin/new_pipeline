
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

    bam_location = datadir / "BAM" / sample / "BAM" / sample
    source_bam = bam_location.with_suffix(f".{bam_file}.bam")
    target_bam = bam_location.with_suffix(".bam")
    source_bai = bam_location.with_suffix(f".{bam_file}.bam.bai")
    target_bai = bam_location.with_suffix(".bam.bai")

    if source_bam.exists():
        console.log(f"BAM file {bam_file}.bam exists")
        console.log(f"Moving {source_bam} to {target_bam}")
        try:
            move(str(source_bam), str(target_bam))
            if source_bai.exists():
                move(str(source_bai), str(target_bai))
            else:
                console.log(f"Index file {source_bai} does not exist. Skipping index move.")
                log_to_api(f"Index file {source_bai} does not exist. Skipping index move.", "WARNING", "move_bam", sample, str(datadir))
            console.log("Files moved")
            log_to_api("Files moved", "INFO", "move_bam", sample, str(datadir))
            return True
        except Exception as e:
            console.log(f"Error moving files: {str(e)}")
            log_to_api(f"Error moving files: {str(e)}", "ERROR", "move_bam", sample, str(datadir))
            return False
    else:
        console.log(f"BAM file {bam_file}.bam does not exist, but process might have been completed")
        log_to_api(
            f"BAM file {bam_file}.bam does not exist, but process might have been completed",
            "WARNING",
            "move_bam",
            sample,
            str(datadir)
        )
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
