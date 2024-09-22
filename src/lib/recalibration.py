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


def recalibrate(datadir: Path, sample_id: str, reference: Path, gatk: str, samtools: str, db: TinyDB) -> str:
    @timer_with_db_log(db)
    def _recalibrate():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        input_bam = bam_dir / f"{sample_id}.bam"
        final_bam = bam_dir / f"{sample_id}.bam"  # This will be our final output
        recal_table = bam_dir / "recal_data.table"
        recalibration_log = bam_dir / "recalibration.txt"

        if not recal_table.exists():
            console.print(f"[bold red]Error: {recal_table} does not exist. Run base_recal1 first.[/bold red]")
            log_to_api("recal_data.table does not exist", "ERROR", "recalibrate", sample_id, str(datadir))
            return "error"

        if not input_bam.exists():
            console.print(f"[bold red]Error: Input BAM file {input_bam} does not exist.[/bold red]")
            log_to_api("Input BAM file does not exist", "ERROR", "recalibrate", sample_id, str(datadir))
            return "error"

        console.print(f"[bold cyan]Starting recalibration for {sample_id}[/bold cyan]")

        # Apply BQSR and output directly to the final BAM name
        GATK_string = (
            f"{gatk} ApplyBQSR -R {reference} -I {input_bam} "
            f"--bqsr-recal-file {recal_table} -O {final_bam}"
        )

        result = run_command(GATK_string, f"Recalibration for {sample_id}", sample_id, datadir)

        if result == 0:
            # Index the recalibrated BAM file
            index_command = f"{samtools} index {final_bam}"
            index_result = run_command(index_command, f"Indexing Recalibrated BAM for {sample_id}", sample_id, datadir)

            if index_result == 0:
                with open(recalibration_log, "w") as recal_file:
                    recal_file.write(f"{GATK_string}\n")
                    recal_file.write(f"{index_command}\n")
                console.print(f"[bold green]Recalibration and indexing completed successfully for {sample_id}[/bold green]")
                log_to_api("Recalibration and indexing completed successfully", "INFO", "recalibrate", sample_id, str(datadir))
                return "success"
            else:
                console.print(f"[bold red]Error: Failed to index recalibrated BAM for {sample_id}[/bold red]")
                log_to_api("Failed to index recalibrated BAM", "ERROR", "recalibrate", sample_id, str(datadir))
                return "error"
        else:
            console.print(f"[bold red]Error: Recalibration failed for {sample_id}[/bold red]")
            log_to_api("Recalibration failed", "ERROR", "recalibrate", sample_id, str(datadir))
            return "error"

    return _recalibrate()
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