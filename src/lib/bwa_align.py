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
def align_reads(sample_id: str, bwa: str, reference: Path, fastq_files: List[Path], threads: int, samtools: str, unsorted_bam: Path, datadir: Path) -> int:
    bam_header = f"@RG\\tID:{sample_id}\\tLB:{datadir}\\tPL:Illumina\\tSM:{sample_id}\\tPU:None"
    bwa_string = (
        f"{bwa} mem -t {threads} -Y -R '{bam_header}' {reference} "
        f"{' '.join([str(file) for file in fastq_files])} | "
        f"{samtools} view -@ {threads // 2} -bS -F 0 - > {unsorted_bam}"
    )
    
    console.print(Panel("[bold cyan]Aligning reads[/bold cyan]"))
    result = run_command(bwa_string, f"Aligning reads for {sample_id}")
    
    if result.returncode != 0:
        console.print(f"[bold red]BWA alignment failed for sample {sample_id}[/bold red]")
        log_to_api(f"BWA alignment failed for sample {sample_id}", "ERROR", "bwa_align", sample_id, str(datadir))
        return 1
    return 0

@timer
def sort_bam(sample_id: str, samtools: str, threads: int, sort_memory: str, temp_dir: Path, unsorted_bam: Path, bam_output: Path) -> int:
    sort_string = (
        f"{samtools} sort -@ {threads // 2} -m {sort_memory} "
        f"-T {temp_dir}/sort_{sample_id} -O bam -o {bam_output} {unsorted_bam}"
    )

    console.print(Panel("[bold cyan]Sorting BAM file[/bold cyan]"))
    result = run_command(sort_string, f"Sorting BAM file for {sample_id}")

    if result.returncode != 0:
        console.print(f"[bold red]BAM sorting failed for sample {sample_id}[/bold red]")
        log_to_api(f"BAM sorting failed for sample {sample_id}", "ERROR", "bwa_align", sample_id, str(datadir))
        return 1
    return 0

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
        sort_memory: str = "4G",
        temp_dir: Path = Path("/tmp")
) -> int:
    console.print(Panel(f"[bold blue]Starting BWA processing for sample: {sample_id}[/bold blue]"))
    log_to_api(
        f"Starting bwa processing for file {sample_id}",
        "INFO",
        "bwa_align",
        sample_id,
        str(datadir),
    )

    if check_existing_bam(sample_id, datadir, threads, samtools) == 0:
        return 0

    console.print(f"[yellow]No existing BAM file found for {sample_id}. Proceeding with alignment.[/yellow]")
    log_to_api(
        f"No existing BAM file found for {sample_id}. Proceeding with alignment.",
        "INFO",
        "bwa_align",
        sample_id,
        str(datadir),
    )

    bam_output = datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.dedup.bam"
    metrics_file = datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.markdup_metrics.txt"

    bam_header = f"@RG\\tID:{sample_id}\\tLB:{datadir}\\tPL:Illumina\\tSM:{sample_id}\\tPU:None"

    # Combine alignment, sorting, duplicate marking, and indexing in one command
    command = (
        f"{bwa} mem -t {threads} -Y -R '{bam_header}' {reference} "
        f"{' '.join([str(file) for file in fastq_files])} | "
        f"{samtools} view -@ {threads // 4} -bS -F 0 - | "
        f"{samtools} sort -@ {threads // 4} -m {sort_memory} -T {temp_dir}/sort_{sample_id} - | "
        f"{samtools} markdup -@ {threads // 4} -s - {bam_output} 2> {metrics_file} && "
        f"{samtools} index -@ {threads} {bam_output}"
    )

    console.print(Panel("[bold cyan]Aligning, sorting, marking duplicates, and indexing[/bold cyan]"))
    result = run_command(command, f"Processing {sample_id}")

    if result.returncode != 0:
        console.print(f"[bold red]Processing failed for sample {sample_id}[/bold red]")
        log_to_api(f"Processing failed for sample {sample_id}", "ERROR", "bwa_align", sample_id, str(datadir))
        return 1

    if bam_output.exists() and bam_output.with_suffix('.bam.bai').exists():
        console.print(f"[bold green]Deduplicated BAM and BAI files created successfully for {sample_id}[/bold green]")
        log_to_api(
            f"Deduplicated BAM and BAI files created successfully for {sample_id}",
            "INFO",
            "bwa_align",
            sample_id,
            str(datadir),
        )
        return 0
    else:
        console.print(f"[bold red]Failed to create deduplicated BAM or BAI file for {sample_id}[/bold red]")
        log_to_api(f"Failed to create deduplicated BAM or BAI file for {sample_id}", "ERROR", "bwa_align", sample_id,
                   str(datadir))
        return 1

def print_summary(timings: dict):
    table = Table(title="BWA Alignment Process Summary")
    table.add_column("Step", style="cyan", no_wrap=True)
    table.add_column("Time (seconds)", style="magenta")

    for step, time in timings.items():
        table.add_row(step, f"{time:.2f}")

    console.print(table)

if __name__ == "__main__":
    pass
