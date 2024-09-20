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
import time
from functools import wraps

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn
from rich.table import Table

from .utils import move_bam
from .log_api import log_to_api

console = Console()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        duration = end_time - start_time
        function_name = func.__name__
        console.print(f"[bold cyan]{function_name} completed in {duration:.2f} seconds[/bold cyan]")
        return result
    return wrapper


@timer
def run_command(command: str, description: str) -> subprocess.CompletedProcess:
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console,
    ) as progress:
        task = progress.add_task(description, total=None)
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        
        while True:
            output = process.stderr.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                progress.console.print(output.strip(), style="dim")
            progress.update(task, advance=1)

        rc = process.poll()
        return subprocess.CompletedProcess(process.args, rc, stdout=process.stdout.read(), stderr=process.stderr.read())

@timer
def check_existing_bam(sample_id: str, datadir: Path, threads: int, samtools: str) -> int:
    bam_paths = [
        datadir / "BAM" / sample_id / f"{sample_id}.bam",
        datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.bam",
        datadir / "BAM" / sample_id / f"{sample_id}.bam.gz",
        datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.bam.gz"
    ]
    
    existing_bam_path = next((path for path in bam_paths if path.exists()), None)

    if existing_bam_path:
        console.print(f"[green]BAM file exists:[/green] {existing_bam_path}")
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
            console.print(f"[yellow]BAI file not found for {sample_id}. Creating index.[/yellow]")
            log_to_api(
                f"BAI file not found for {sample_id}. Creating index.",
                "INFO",
                "bwa_align",
                sample_id,
                str(datadir),
            )
            index_string = f"{samtools} index -@ {threads} {existing_bam_path}"
            result = run_command(index_string, f"Indexing BAM file for {sample_id}")
            
            if bai_path.exists():
                console.print(f"[green]BAI file created for {sample_id}[/green]")
                log_to_api(
                    f"BAI file created for {sample_id}",
                    "INFO",
                    "bwa_align",
                    sample_id,
                    str(datadir),
                )
            else:
                console.print(f"[bold red]Failed to create BAI file for {sample_id}[/bold red]")
                log_to_api(
                    f"Failed to create BAI file for {sample_id}",
                    "ERROR",
                    "bwa_align",
                    sample_id,
                    str(datadir),
                )
        return 0
    return 1


@timer
def index_bam(sample_id: str, samtools: str, threads: int, bam_output: Path) -> int:
    index_string = f"{samtools} index -@ {threads} {bam_output}"
    console.print(Panel("[bold cyan]Indexing BAM file[/bold cyan]"))
    result = run_command(index_string, f"Indexing BAM file for {sample_id}")

    if result.returncode != 0:
        console.print(f"[bold red]BAM indexing failed for sample {sample_id}[/bold red]")
        log_to_api(f"BAM indexing failed for sample {sample_id}", "ERROR", "bwa_align", sample_id, str(datadir))
        return 1
    return 0


@timer
def run_bwa(
        sample_id: str,
        fastq_files: List[Path],
        datadir: Path,
        reference: Path,
        bwa: str,
        samtools: str,
        threads: int = 16,
        sort_memory: str = "4G"
) -> int:
    bam_dir = datadir / "BAM" / sample_id / "BAM"
    bam_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = bam_dir / "tmp"
    temp_dir.mkdir(exist_ok=True)

    bam_output = bam_dir / f"{sample_id}.bam"
    metrics_file = bam_dir / f"{sample_id}.markdup_metrics.txt"

    # Check if BAM file already exists
    if bam_output.exists():
        file_size = bam_output.stat().st_size
        console.print(f"[yellow]BAM file already exists for {sample_id}. Size: {file_size} bytes[/yellow]")
        if file_size > 1_000_000:  # More than 1MB
            console.print("[green]Existing BAM file seems valid. Skipping alignment process.[/green]")
            return 0
        else:
            console.print("[yellow]Existing BAM file is suspiciously small. Will recreate.[/yellow]")

    bam_header = f"@RG\\tID:{sample_id}\\tLB:{datadir}\\tPL:Illumina\\tSM:{sample_id}\\tPU:None"

    steps = [
        # Alignment
        (f"{bwa} mem -t {threads} -Y -R '{bam_header}' {reference} "
         f"{' '.join([str(file) for file in fastq_files])} > {temp_dir}/{sample_id}.sam",
         f"BWA alignment for {sample_id}"),

        # Convert SAM to BAM
        (f"{samtools} view -@ {threads // 2} -bh {temp_dir}/{sample_id}.sam > {temp_dir}/{sample_id}.bam",
         f"Convert SAM to BAM for {sample_id}"),

        # Sort BAM by name (required for fixmate)
        (f"{samtools} sort -@ {threads // 2} -n -o {temp_dir}/{sample_id}.namesort.bam {temp_dir}/{sample_id}.bam",
         f"Sort BAM by name for {sample_id}"),

        # Fixmate
        (
        f"{samtools} fixmate -@ {threads // 2} -m {temp_dir}/{sample_id}.namesort.bam {temp_dir}/{sample_id}.fixmate.bam",
        f"Fixmate for {sample_id}"),

        # Sort BAM by position
        (f"{samtools} sort -@ {threads // 2} -m {sort_memory} -T {temp_dir}/sort_{sample_id} "
         f"-o {temp_dir}/{sample_id}.sorted.bam {temp_dir}/{sample_id}.fixmate.bam",
         f"Sort BAM by position for {sample_id}"),

        # Mark duplicates
        (f"{samtools} markdup -@ {threads // 2} -s {temp_dir}/{sample_id}.sorted.bam {bam_output} "
         f"2> {metrics_file}",
         f"Mark duplicates for {sample_id}"),

        # Index BAM
        (f"{samtools} index -@ {threads} {bam_output}",
         f"Index BAM for {sample_id}")
    ]

    for command, description in steps:
        console.print(f"[cyan]Executing: {description}[/cyan]")
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            console.print(f"[bold red]Error in {description}[/bold red]")
            console.print(f"Command: {command}")
            console.print(f"STDOUT: {result.stdout}")
            console.print(f"STDERR: {result.stderr}")
            return 1

        # Check file sizes after each step
        output_file = command.split(">")[-1].strip().split()[0] if ">" in command else command.split()[-1]
        if output_file.endswith('.bam') or output_file.endswith('.sam'):
            file_size = Path(output_file).stat().st_size
            console.print(f"[green]File size of {output_file}: {file_size} bytes[/green]")

    # Clean up temporary files
    shutil.rmtree(temp_dir)

    # Final file size check
    final_size = bam_output.stat().st_size
    console.print(f"[bold green]Final BAM file size: {final_size} bytes[/bold green]")

    if final_size < 1_000_000:  # Less than 1MB
        console.print("[bold red]Warning: Final BAM file is suspiciously small![/bold red]")
        return 1

    return 0


def print_summary(timings: dict):
    table = Table(title="BWA Alignment Process Summary")
    table.add_column("Step", style="cyan", no_wrap=True)
    table.add_column("Time (seconds)", style="magenta")

    for step, time in timings.items():
        table.add_row(step, f"{time:.2f}")

    console.print(table)


if __name__ == "__main__":
    pass
