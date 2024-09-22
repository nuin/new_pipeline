# lib/recalibration.py

import subprocess
from pathlib import Path
from typing import Dict
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax
from tinydb import TinyDB
from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log

console = Console()


def run_command(command: str, description: str, sample_id: str, datadir: Path) -> int:
    console.print(Panel(f"[bold blue]{description}[/bold blue]"))
    console.print(Syntax(command, "bash", theme="monokai", line_numbers=True))

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                               bufsize=1, universal_newlines=True)

    output = []
    with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console,
            transient=True
    ) as progress:
        task = progress.add_task(description, total=None)

        for line in process.stdout:
            output.append(line)
            progress.update(task, advance=1)
            console.print(f"[dim]{line.strip()}[/dim]")

    process.wait()

    if process.returncode != 0:
        console.print(f"[bold red]Error in {description}[/bold red]")
        console.print("Command output:")
        console.print("".join(output))
        log_to_api(f"Error in {description}", "ERROR", "recalibration", sample_id, str(datadir))
        return 1

    log_to_api(f"{description} completed successfully", "INFO", "recalibration", sample_id, str(datadir))
    return 0


def base_recal1(datadir: Path, sample_id: str, bed_file: Path, vcf_file: Path, reference: Path, gatk: str,
                db: TinyDB) -> str:
    @timer_with_db_log(db)
    def _base_recal1():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        recal_table = bam_dir / "recal_data.table"

        if recal_table.exists():
            console.print(f"[yellow]{recal_table} file exists. Skipping base recalibration step 1.[/yellow]")
            log_to_api("recal_data.table file exists", "INFO", "base_recal1", sample_id, str(datadir))
            return "exists"

        console.print(f"[bold cyan]Starting base recalibration step 1 for {sample_id}[/bold cyan]")

        GATK_string = (
            f"{gatk} BaseRecalibrator -R {reference} "
            f"-I {bam_dir}/{sample_id}.bam --known-sites {vcf_file} "
            f"-O {recal_table} -L {bed_file}"
        )

        result = run_command(GATK_string, f"Base Recalibration Step 1 for {sample_id}", sample_id, datadir)

        return "success" if result == 0 else "error"

    return _base_recal1()


import tempfile
import shutil

import tempfile
import shutil
import subprocess
from pathlib import Path


def _recalibrate():
    bam_dir = datadir / "BAM" / sample_id / "BAM"
    input_bam = bam_dir / f"{sample_id}.bam"
    recal_table = bam_dir / "recal_data.table"
    recalibration_log = bam_dir / "recalibration.log"

    # Check if input files exist
    if not input_bam.exists():
        with open(recalibration_log, 'w') as log_file:
            log_file.write(f"Error: Input BAM file not found: {input_bam}\n")
        return False

    if not recal_table.exists():
        with open(recalibration_log, 'w') as log_file:
            log_file.write(f"Error: Recalibration table not found: {recal_table}\n")
        return False

    # Create a temporary directory for the output
    with tempfile.TemporaryDirectory(dir=bam_dir) as temp_dir:
        temp_output_bam = Path(temp_dir) / f"{sample_id}.recal.bam"

        # Construct the GATK command
        command = f"{gatk} ApplyBQSR " \
                  f"-R {reference} " \
                  f"-I {input_bam} " \
                  f"--bqsr-recal-file {recal_table} " \
                  f"-O {temp_output_bam} " \
                  f"--create-output-bam-index true"

        # Run the recalibration
        try:
            result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)

            # Write the log
            with open(recalibration_log, 'w') as log_file:
                log_file.write(f"Command executed: {command}\n\n")
                log_file.write("STDOUT:\n")
                log_file.write(result.stdout)
                log_file.write("\nSTDERR:\n")
                log_file.write(result.stderr)

            # Check if the output BAM was created
            if not temp_output_bam.exists():
                with open(recalibration_log, 'a') as log_file:
                    log_file.write(f"\nError: Output BAM file was not created: {temp_output_bam}\n")
                return False

            # If successful, replace the original BAM file
            shutil.move(str(temp_output_bam), str(input_bam))

            # Also move the index file if it was created
            temp_bai = temp_output_bam.with_suffix('.bai')
            if temp_bai.exists():
                shutil.move(str(temp_bai), str(input_bam.with_suffix('.bai')))

            with open(recalibration_log, 'a') as log_file:
                log_file.write("\nRecalibration completed successfully.\n")

            return True

        except subprocess.CalledProcessError as e:
            # Log the error
            with open(recalibration_log, 'w') as log_file:
                log_file.write(f"Command failed: {command}\n")
                log_file.write(f"Error: {str(e)}\n")
                log_file.write("STDOUT:\n")
                log_file.write(e.stdout)
                log_file.write("\nSTDERR:\n")
                log_file.write(e.stderr)
            return False

        except Exception as e:
            # Log any other unexpected errors
            with open(recalibration_log, 'w') as log_file:
                log_file.write(f"Unexpected error occurred: {str(e)}\n")
            return False

    return False  # If we get here, something went wrong


def is_bam_recalibrated(bam_path: Path) -> bool:
    """
    Check if the BAM file has already been recalibrated.
    """
    recal_bam = bam_path.with_name(f"{bam_path.stem}.recal_reads.bam")
    recal_log = bam_path.with_name("recalibration.txt")
    return recal_bam.exists() and recal_log.exists()


def recalibration_pipeline(datadir: Path, sample_id: str, bed_file: Path, vcf_file: Path, reference: Path, gatk: str,
                           samtools: str, db: TinyDB) -> Dict[str, int]:
    @timer_with_db_log(db)
    def _recalibration_pipeline():
        log_to_db(db, "Starting recalibration pipeline", "INFO", "recalibration", sample_id, datadir.name)

        bam_dir = datadir / "BAM" / sample_id / "BAM"
        input_bam = bam_dir / f"{sample_id}.bam"

        if is_bam_recalibrated(input_bam):
            console.print(f"[yellow]BAM file for {sample_id} is already recalibrated. Skipping recalibration.[/yellow]")
            log_to_api(f"BAM file for {sample_id} is already recalibrated", "INFO", "recalibration", sample_id, str(datadir))
            log_to_db(db, f"BAM file for {sample_id} is already recalibrated", "INFO", "recalibration", sample_id, datadir.name)
            return {"status": 0}

        recal1_result = base_recal1(datadir, sample_id, bed_file, vcf_file, reference, gatk, db)
        if recal1_result != "success":
            log_to_db(db, "Base recalibration step 1 failed", "ERROR", "recalibration", sample_id, datadir.name)
            return {"status": 1}

        recal2_result = recalibrate(datadir, sample_id, reference, gatk, samtools, db)
        if recal2_result != "success":
            log_to_db(db, "Recalibration step 2 failed", "ERROR", "recalibration", sample_id, datadir.name)
            return {"status": 1}

        log_to_db(db, "Recalibration pipeline completed successfully", "INFO", "recalibration", sample_id, datadir.name)
        return {"status": 0}

    return _recalibration_pipeline()


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass