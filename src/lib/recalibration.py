import subprocess
from pathlib import Path
from typing import Dict
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax
from tinydb import TinyDB

from .log_api import log_to_api
from .db_logger import get_sample_db, log_to_db, timer_with_db_log

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
        log_to_api(f"Error in {description}", "ERROR", "recalibration", sample_id, datadir.name)
        return 1

    log_to_api(f"{description} completed successfully", "INFO", "recalibration", sample_id, datadir.name)
    return 0

@timer_with_db_log(get_sample_db(datadir, sample_id))
def base_recal1(datadir: Path, sample_id: str, bed_file: Path, vcf_file: Path, reference: Path, gatk: str) -> Dict[str, int]:
    bam_dir = datadir / "BAM" / sample_id / "BAM"
    recal_table = bam_dir / "recal_data.table"

    if recal_table.exists():
        console.print(f"[yellow]{recal_table} file exists. Skipping base recalibration step 1.[/yellow]")
        log_to_api("recal_data.table file exists", "INFO", "base_recal1", sample_id, datadir.name)
        return {"status": 0}

    console.print(f"[bold cyan]Starting base recalibration step 1 for {sample_id}[/bold cyan]")

    GATK_string = (
        f"{gatk} BaseRecalibrator -R {reference} "
        f"-I {bam_dir}/{sample_id}.bam --known-sites {vcf_file} "
        f"-O {recal_table} -L {bed_file}"
    )

    result = run_command(GATK_string, f"Base Recalibration Step 1 for {sample_id}", sample_id, datadir)

    return {"status": result}

@timer_with_db_log(get_sample_db(datadir, sample_id))
def recalibrate(datadir: Path, sample_id: str, reference: Path, gatk: str, samtools: str) -> Dict[str, int]:
    bam_dir = datadir / "BAM" / sample_id / "BAM"
    recal_bam = bam_dir / f"{sample_id}.recal_reads.bam"
    recal_table = bam_dir / "recal_data.table"
    recalibration_log = bam_dir / "recalibration.txt"

    if recal_bam.exists():
        console.print(f"[yellow]{recal_bam} file exists. Skipping recalibration step 2.[/yellow]")
        log_to_api("recal_reads.bam file exists", "INFO", "recalibrate", sample_id, datadir.name)
        return {"status": 0}
    
    if not recal_table.exists():
        console.print(f"[bold red]Error: {recal_table} does not exist. Run base_recal1 first.[/bold red]")
        log_to_api("recal_data.table does not exist", "ERROR", "recalibrate", sample_id, datadir.name)
        return {"status": 1}

    console.print(f"[bold cyan]Starting recalibration step 2 for {sample_id}[/bold cyan]")

    GATK_string = (
        f"{gatk} ApplyBQSR -R {reference} -I {bam_dir}/{sample_id}.bam "
        f"--bqsr-recal-file {recal_table} -O {recal_bam}"
    )

    result = run_command(GATK_string, f"Recalibration Step 2 for {sample_id}", sample_id, datadir)

    if result == 0:
        # Index the recalibrated BAM file
        index_command = f"{samtools} index {recal_bam}"
        index_result = run_command(index_command, f"Indexing Recalibrated BAM for {sample_id}", sample_id, datadir)

        if index_result == 0:
            with open(recalibration_log, "w") as recal_file:
                recal_file.write(f"{GATK_string}\n")
                recal_file.write(f"{index_command}\n")
        else:
            result = 1  # Set result to error if indexing fails

    return {"status": result}

def recalibration_pipeline(datadir: Path, sample_id: str, bed_file: Path, vcf_file: Path, reference: Path, gatk: str, samtools: str) -> Dict[str, int]:
    db = get_sample_db(datadir, sample_id)

    log_to_db(db, "Starting recalibration pipeline", "INFO", "recalibration", sample_id, datadir.name)

    recal1_result = base_recal1(datadir, sample_id, bed_file, vcf_file, reference, gatk)
    if recal1_result["status"] != 0:
        log_to_db(db, "Base recalibration step 1 failed", "ERROR", "recalibration", sample_id, datadir.name)
        return {"status": 1}

    recal2_result = recalibrate(datadir, sample_id, reference, gatk, samtools)
    if recal2_result["status"] != 0:
        log_to_db(db, "Recalibration step 2 failed", "ERROR", "recalibration", sample_id, datadir.name)
        return {"status": 1}

    log_to_db(db, "Recalibration pipeline completed successfully", "INFO", "recalibration", sample_id, datadir.name)
    return {"status": 0}

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass