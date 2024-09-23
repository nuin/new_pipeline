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
import time
import traceback

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
            recal_table = bam_dir / "recal_data.table"

            # Log the existence of files
            console.print(f"[cyan]Checking for existing files:[/cyan]")
            console.print(f"Input BAM: {input_bam.exists()}")
            console.print(f"BAM index: {input_bam.with_suffix('.bam.bai').exists()}")
            console.print(f"Recal table: {recal_table.exists()}")
            console.print(f"Recalibration flag: {recalibration_flag.exists()}")

            # Check if recalibration was already completed
            if input_bam.exists() and input_bam.with_suffix('.bam.bai').exists() and recal_table.exists():
                console.print(
                    f"[yellow]Recalibrated BAM and index exist for {sample_id}. Assuming recalibration was completed.[/yellow]")
                log_to_api(f"Recalibrated BAM and index exist for {sample_id}", "INFO", "recalibration", sample_id,
                           str(datadir))
                log_to_db(db, f"Recalibrated BAM and index exist for {sample_id}", "INFO", "recalibration", sample_id,
                          datadir.name)

                # Create the flag file if it doesn't exist
                if not recalibration_flag.exists():
                    try:
                        recalibration_flag.touch()
                        console.print(f"[green]Recalibration flag created for {sample_id}[/green]")
                        log_to_db(db, f"Recalibration flag created for {sample_id}", "INFO", "recalibration", sample_id,
                                  datadir.name)
                    except Exception as e:
                        console.print(
                            f"[bold yellow]Warning: Failed to create recalibration flag: {str(e)}[/bold yellow]")
                        log_to_db(db, f"Warning: Failed to create recalibration flag: {str(e)}", "WARNING",
                                  "recalibration", sample_id, datadir.name)

                # Even if flag creation fails, we assume recalibration is complete
                log_to_db(db, "Recalibration pipeline completed successfully (pre-existing files)", "INFO",
                          "recalibration", sample_id, datadir.name)
                console.print(
                    f"[bold green]Recalibration pipeline completed successfully for {sample_id} (pre-existing files)[/bold green]")
                return {"status": 0}

            # If we reach here, recalibration needs to be performed
            if not recal_table.exists():
                console.print(f"[cyan]Starting base recalibration step 1 (BaseRecalibrator) for {sample_id}[/cyan]")
                recal1_result = base_recal1(datadir, sample_id, bed_file, vcf_file, reference, gatk, db)
                console.print(f"[cyan]Base recalibration step 1 result: {recal1_result}[/cyan]")

                if recal1_result != "success":
                    error_msg = f"Base recalibration step 1 failed for {sample_id}: {recal1_result}"
                    log_to_db(db, error_msg, "ERROR", "recalibration", sample_id, datadir.name)
                    console.print(f"[bold red]{error_msg}[/bold red]")
                    return {"status": 1}
            else:
                console.print(
                    f"[yellow]Recalibration table exists for {sample_id}. Skipping BaseRecalibrator step.[/yellow]")

            console.print(f"[cyan]Starting recalibration step 2 (ApplyBQSR) for {sample_id}[/cyan]")
            log_to_db(db, "Starting recalibration step 2 (ApplyBQSR)", "INFO", "recalibration", sample_id, datadir.name)

            recal2_result = recalibrate(datadir, sample_id, reference, gatk, samtools, db)
            console.print(f"[cyan]Recalibration step 2 (ApplyBQSR) result: {recal2_result}[/cyan]")

            if recal2_result != "success":
                error_msg = f"Recalibration step 2 (ApplyBQSR) failed for {sample_id}: {recal2_result}"
                log_to_db(db, error_msg, "ERROR", "recalibration", sample_id, datadir.name)
                console.print(f"[bold red]{error_msg}[/bold red]")
                return {"status": 1}

            # Create the recalibration flag file
            try:
                recalibration_flag.touch()
                console.print(f"[green]Recalibration flag created for {sample_id}[/green]")
                log_to_db(db, f"Recalibration flag created for {sample_id}", "INFO", "recalibration", sample_id,
                          datadir.name)
            except Exception as e:
                console.print(f"[bold yellow]Warning: Failed to create recalibration flag: {str(e)}[/bold yellow]")
                log_to_db(db, f"Warning: Failed to create recalibration flag: {str(e)}", "WARNING", "recalibration",
                          sample_id, datadir.name)

            log_to_db(db, "Recalibration pipeline completed successfully", "INFO", "recalibration", sample_id,
                      datadir.name)
            console.print(f"[bold green]Recalibration pipeline completed successfully for {sample_id}[/bold green]")
            return {"status": 0}
        except Exception as e:
            error_msg = f"Unexpected error in recalibration pipeline for {sample_id}: {str(e)}\n"
            error_msg += f"Traceback:\n{''.join(traceback.format_tb(e.__traceback__))}"
            log_to_db(db, error_msg, "ERROR", "recalibration", sample_id, datadir.name)
            console.print(f"[bold red]{error_msg}[/bold red]")
            return {"status": 1}

    return _recalibration_pipeline()


def recalibrate(datadir: Path, sample_id: str, reference: Path, gatk: str, samtools: str, db: TinyDB) -> str:
    @timer_with_db_log(db)
    def _recalibrate():
        try:
            bam_dir = datadir / "BAM" / sample_id / "BAM"
            input_bam = bam_dir / f"{sample_id}.bam"
            recal_table = bam_dir / "recal_data.table"
            output_bam = bam_dir / f"{sample_id}.recalibrated.bam"

            console.print(Panel(f"[bold blue]Starting ApplyBQSR for {sample_id}[/bold blue]"))
            log_to_db(db, f"Starting ApplyBQSR for {sample_id}", "INFO", "recalibration", sample_id, datadir.name)

            command = f"{gatk} ApplyBQSR " \
                      f"-R {reference} " \
                      f"-I {input_bam} " \
                      f"--bqsr-recal-file {recal_table} " \
                      f"-O {output_bam}"

            console.print(Panel(f"[cyan]Executing command:[/cyan]\n{command}"))
            log_to_db(db, f"Executing ApplyBQSR command: {command}", "INFO", "recalibration", sample_id, datadir.name)

            with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    TimeElapsedColumn(),
                    console=console
            ) as progress:
                task = progress.add_task("[cyan]Running ApplyBQSR...", total=None)

                process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                           universal_newlines=True)

                output = []
                for line in process.stdout:
                    output.append(line)
                    progress.update(task, advance=1)
                    log_to_db(db, line.strip(), "DEBUG", "recalibration", sample_id, datadir.name)

                process.wait()

                if process.returncode != 0:
                    raise subprocess.CalledProcessError(process.returncode, command, output=''.join(output))

            if not output_bam.exists():
                raise FileNotFoundError(f"Output BAM file was not created: {output_bam}")

            console.print(Panel(f"[green]ApplyBQSR completed successfully for {sample_id}[/green]"))
            log_to_db(db, f"ApplyBQSR completed successfully for {sample_id}", "INFO", "recalibration", sample_id,
                      datadir.name)

            # Replace the original BAM with the recalibrated one
            console.print(f"[cyan]Replacing original BAM with recalibrated BAM for {sample_id}[/cyan]")
            input_bam.unlink()
            output_bam.rename(input_bam)

            time.sleep(2)  # Wait for 2 seconds to ensure file system operations are complete

            # Create BAM index
            console.print(f"[cyan]Creating BAM index for {sample_id}[/cyan]")
            index_command = f"{samtools} index {input_bam}"

            with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    TimeElapsedColumn(),
                    console=console
            ) as progress:
                task = progress.add_task("[cyan]Indexing BAM file...", total=None)

                index_process = subprocess.Popen(index_command, shell=True, stdout=subprocess.PIPE,
                                                 stderr=subprocess.STDOUT, universal_newlines=True)

                index_output = []
                while True:
                    output = index_process.stdout.readline()
                    if output == '' and index_process.poll() is not None:
                        break
                    if output:
                        index_output.append(output)
                        progress.update(task, advance=1)
                        log_to_db(db, output.strip(), "DEBUG", "recalibration", sample_id, datadir.name)

                index_process.wait()

                if index_process.returncode != 0:
                    raise subprocess.CalledProcessError(index_process.returncode, index_command,
                                                        output=''.join(index_output))

            bai_file = input_bam.with_suffix('.bam.bai')
            if not bai_file.exists():
                raise FileNotFoundError(f"BAM index file was not created: {bai_file}")

            console.print(
                Panel(f"[bold green]ApplyBQSR and indexing completed successfully for {sample_id}[/bold green]"))
            log_to_db(db, f"ApplyBQSR and indexing completed successfully for {sample_id}", "INFO", "recalibration",
                      sample_id, datadir.name)

            return "success"

        except subprocess.CalledProcessError as e:
            error_msg = f"Error in ApplyBQSR or indexing for {sample_id}:\n"
            error_msg += f"Command: {e.cmd}\n"
            error_msg += f"Return code: {e.returncode}\n"
            error_msg += f"Output: {e.output}\n"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_db(db, error_msg, "ERROR", "recalibration", sample_id, datadir.name)
            return f"error: {error_msg}"

        except FileNotFoundError as e:
            error_msg = f"File not found error in recalibration for {sample_id}: {e}\n"
            error_msg += f"Current working directory: {os.getcwd()}\n"
            error_msg += f"Directory contents: {os.listdir(bam_dir)}\n"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_db(db, error_msg, "ERROR", "recalibration", sample_id, datadir.name)
            return f"error: {error_msg}"

        except Exception as e:
            error_msg = f"Unexpected error in ApplyBQSR or indexing for {sample_id}: {e}\n"
            error_msg += f"Traceback:\n{''.join(traceback.format_tb(e.__traceback__))}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_db(db, error_msg, "ERROR", "recalibration", sample_id, datadir.name)
            return f"error: {error_msg}"

    return _recalibrate()

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass