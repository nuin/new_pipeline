"""
.. module:: bwa_align
    :platform: Any
    :synopsis: Module that performs BWA alignment and generates BAM files
.. moduleauthor:: Paulo Nuin, July 2016
"""

import logging
import subprocess
import shutil
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
    fastq_files: List[Path],
    datadir: Path,
    reference: Path,
    bwa: str,
    samtools: str,
    threads: int = 16,
    sort_memory: str = "4G",
    temp_dir: Path = Path("/tmp")
) -> int:
    """
    Performs the alignment/mapping for each pair/sample_id.

    Parameters:
    sample_id (str): The ID of the sample to be analysed.
    fastq_files (List[Path]): The list of FASTQ files for the sample.
    datadir (Path): The directory of the run.
    reference (Path): The reference genome.
    bwa (str): The path to the BWA software.
    samtools (str): The path to the Samtools software.
    threads (int): Number of threads to use for alignment and sorting.
    sort_memory (str): Amount of memory to use for sorting.
    temp_dir (Path): Directory for temporary files.

    Returns:
    int: A flag indicating whether the BAM file exists (0) or was created (1).
    """

    console.log(f"Starting bwa processing for file {sample_id}")
    log_to_api(
        f"Starting bwa processing for file {sample_id}",
        "INFO",
        "bwa_align",
        sample_id,
        str(datadir),
    )

    bam_paths = [
        datadir / "BAM" / sample_id / f"{sample_id}.bam",
        datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.bam",
        datadir / "BAM" / sample_id / f"{sample_id}.bam.gz",
        datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.bam.gz"
    ]
    
    existing_bam_path = next((path for path in bam_paths if path.exists()), None)

    if existing_bam_path:
        console.log(f"{sample_id} BAM file exists at {existing_bam_path}")
        log_to_api(
            f"{sample_id} BAM file exists at {existing_bam_path}",
            "INFO",
            "bwa_align",
            sample_id,
            str(datadir),
        )
        
        # Check for BAI file
        bai_path = existing_bam_path.with_suffix('.bam.bai')
        if not bai_path.exists():
            console.log(f"BAI file not found for {sample_id}. Creating index.")
            log_to_api(
                f"BAI file not found for {sample_id}. Creating index.",
                "INFO",
                "bwa_align",
                sample_id,
                str(datadir),
            )
            index_string = f"{samtools} index -@ {threads} {existing_bam_path}"
            subprocess.run(index_string, shell=True, check=True)
            
            if bai_path.exists():
                console.log(f"BAI file created for {sample_id}")
                log_to_api(
                    f"BAI file created for {sample_id}",
                    "INFO",
                    "bwa_align",
                    sample_id,
                    str(datadir),
                )
            else:
                console.log(f"Failed to create BAI file for {sample_id}", style="bold red")
                log_to_api(
                    f"Failed to create BAI file for {sample_id}",
                    "ERROR",
                    "bwa_align",
                    sample_id,
                    str(datadir),
                )
        return 0

    # If we reach this point, no existing BAM file was found
    console.log(f"No existing BAM file found for {sample_id}. Proceeding with alignment.")
    log_to_api(
        f"No existing BAM file found for {sample_id}. Proceeding with alignment.",
        "INFO",
        "bwa_align",
        sample_id,
        str(datadir),
    )

    bam_header = f"@RG\\tID:{sample_id}\\tLB:{datadir}\\tPL:Illumina\\tSM:{sample_id}\\tPU:None"
    bam_output = datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.bam"
    unsorted_bam = temp_dir / f"{sample_id}.unsorted.bam"
    
    # Alignment step
    bwa_string = (
        f"{bwa} mem -t {threads} -Y -R '{bam_header}' {reference} "
        f"{' '.join([str(file) for file in fastq_files])} | "
        f"{samtools} view -@ {threads // 2} -bS -F 0 - > {unsorted_bam}"
    )
    
    console.log(f"Aligning: {bwa_string}")
    log_to_api(
        f"Aligning: {bwa_string}",
        "INFO",
        "bwa_align",
        sample_id,
        str(datadir),
    )
    
    try:
        subprocess.run(bwa_string, shell=True, check=True)
    except subprocess.CalledProcessError:
        error_message = f"BWA alignment failed for sample {sample_id}"
        console.log(error_message, style="bold red")
        log_to_api(error_message, "ERROR", "bwa_align", sample_id, str(datadir))
        return 1

    # Sorting step
    sort_string = (
        f"{samtools} sort -@ {threads // 2} -m {sort_memory} "
        f"-T {temp_dir}/sort_{sample_id} -O bam -o {bam_output} {unsorted_bam}"
    )

    console.log(f"Sorting: {sort_string}")
    log_to_api(
        f"Sorting: {sort_string}",
        "INFO",
        "bwa_align",
        sample_id,
        str(datadir),
    )

    try:
        subprocess.run(sort_string, shell=True, check=True)
    except subprocess.CalledProcessError:
        error_message = f"BAM sorting failed for sample {sample_id}"
        console.log(error_message, style="bold red")
        log_to_api(error_message, "ERROR", "bwa_align", sample_id, str(datadir))
        return 1

    # Remove unsorted BAM
    unsorted_bam.unlink()

    # Create index for the BAM file
    index_string = f"{samtools} index -@ {threads} {bam_output}"
    try:
        subprocess.run(index_string, shell=True, check=True)
    except subprocess.CalledProcessError:
        error_message = f"BAM indexing failed for sample {sample_id}"
        console.log(error_message, style="bold red")
        log_to_api(error_message, "ERROR", "bwa_align", sample_id, str(datadir))
        return 1

    if bam_output.exists() and bam_output.with_suffix('.bam.bai').exists():
        console.log(f"BAM and BAI files created successfully for {sample_id}")
        log_to_api(
            f"BAM and BAI files created successfully for {sample_id}",
            "INFO",
            "bwa_align",
            sample_id,
            str(datadir),
        )
        return 0
    else:
        error_message = f"Failed to create BAM or BAI file for {sample_id}"
        console.log(error_message, style="bold red")
        log_to_api(error_message, "ERROR", "bwa_align", sample_id, str(datadir))
        return 1


if __name__ == "__main__":
    pass
