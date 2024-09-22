# lib/recalibration.py

import subprocess
import tempfile
import shutil
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
        try:
            bam_dir = datadir / "BAM" / sample_id / "BAM"
            input_bam = bam_dir / f"{sample_id}.bam"
            recal_table = bam_dir / "recal_data.table"
            output_bam = bam_dir / f"{sample_id}.recalibrated.bam"

            console.print(f"[cyan]Starting ApplyBQSR for {sample_id}[/cyan]")
            log_to_db(db, f"Starting ApplyBQSR for {sample_id}", "INFO", "recalibration", sample_id, datadir.name)

            command = f"{gatk} ApplyBQSR " \
                      f"-R {reference} " \
                      f"-I {input_bam} " \
                      f"--bqsr-recal-file {recal_table} " \
                      f"-O {output_bam}"

            console.print(f"[dim]Executing command: {command}[/dim]")
            log_to_db(db, f"Executing ApplyBQSR command: {command}", "INFO", "recalibration", sample_id, datadir.name)

            result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)

            if not output_bam.exists():
                raise FileNotFoundError(f"Output BAM file was not created: {output_bam}")

            # Replace the original BAM with the recalibrated one
            input_bam.unlink()
            output_bam.rename(input_bam)

            console.print(f"[green]ApplyBQSR completed successfully for {sample_id}[/green]")
            log_to_db(db, f"ApplyBQSR completed successfully for {sample_id}", "INFO", "recalibration", sample_id, datadir.name)

            return "success"

        except subprocess.CalledProcessError as e:
            console.print(f"[bold red]Error in ApplyBQSR for {sample_id}: {e}[/bold red]")
            log_to_db(db, f"Error in ApplyBQSR: {str(e)}", "ERROR", "recalibration", sample_id, datadir.name)
            return "error"

        except Exception as e:
            console.print(f"[bold red]Unexpected error in ApplyBQSR for {sample_id}: {e}[/bold red]")
            log_to_db(db, f"Unexpected error in ApplyBQSR: {str(e)}", "ERROR", "recalibration", sample_id, datadir.name)
            return "error"

    return _recalibrate()



def is_bam_recalibrated(bam_dir: Path) -> bool:
    recalibration_flag = bam_dir / "recalibration_completed.flag"
    return recalibration_flag.exists()


def recalibration_pipeline(datadir: Path, sample_id: str, bed_file: Path, vcf_file: Path, reference: Path, gatk: str,
                           samtools: str, db: TinyDB) -> Dict[str, int]:
    @timer_with_db_log(db)
    def _recalibration_pipeline():
        try:
            log_to_db(db, "Starting recalibration pipeline", "INFO", "recalibration", sample_id, datadir.name)
            console.print(f"[bold cyan]Starting recalibration pipeline for {sample_id}[/bold cyan]")

            bam_dir = datadir / "BAM" / sample_id / "BAM"
            input_bam = bam_dir / f"{sample_id}.bam"
            recalibration_flag = bam_dir / "recalibration_completed.flag"

            if is_bam_recalibrated(bam_dir):
                console.print(
                    f"[yellow]BAM file for {sample_id} is already recalibrated. Skipping recalibration.[/yellow]")
                log_to_api(f"BAM file for {sample_id} is already recalibrated", "INFO", "recalibration", sample_id,
                           str(datadir))
                log_to_db(db, f"BAM file for {sample_id} is already recalibrated", "INFO", "recalibration", sample_id,
                          datadir.name)
                return {"status": 0}

            console.print(f"[cyan]Starting base recalibration step 1 (BaseRecalibrator) for {sample_id}[/cyan]")
            recal1_result = base_recal1(datadir, sample_id, bed_file, vcf_file, reference, gatk, db)
            console.print(f"[cyan]Base recalibration step 1 result: {recal1_result}[/cyan]")

            if recal1_result != "success":
                log_to_db(db, "Base recalibration step 1 failed", "ERROR", "recalibration", sample_id, datadir.name)
                console.print(f"[bold red]Base recalibration step 1 failed for {sample_id}[/bold red]")
                return {"status": 1}

            console.print(f"[cyan]Starting recalibration step 2 (ApplyBQSR) for {sample_id}[/cyan]")
            log_to_db(db, "Starting recalibration step 2 (ApplyBQSR)", "INFO", "recalibration", sample_id, datadir.name)

            recal2_result = recalibrate(datadir, sample_id, reference, gatk, samtools, db)
            console.print(f"[cyan]Recalibration step 2 (ApplyBQSR) result: {recal2_result}[/cyan]")

            if recal2_result != "success":
                log_to_db(db, "Recalibration step 2 (ApplyBQSR) failed", "ERROR", "recalibration", sample_id,
                          datadir.name)
                console.print(f"[bold red]Recalibration step 2 (ApplyBQSR) failed for {sample_id}[/bold red]")
                return {"status": 1}

            # Create the recalibration flag file
            recalibration_flag.touch()

            log_to_db(db, "Recalibration pipeline completed successfully", "INFO", "recalibration", sample_id,
                      datadir.name)
            console.print(f"[bold green]Recalibration pipeline completed successfully for {sample_id}[/bold green]")
            return {"status": 0}
        except Exception as e:
            log_to_db(db, f"Unexpected error in recalibration pipeline: {str(e)}", "ERROR", "recalibration", sample_id,
                      datadir.name)
            console.print(f"[bold red]Unexpected error in recalibration pipeline for {sample_id}: {str(e)}[/bold red]")
            return {"status": 1}

    return _recalibration_pipeline()


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass