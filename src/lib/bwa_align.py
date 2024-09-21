# lib/bwa_align.py

import logging
import subprocess
import shutil
from pathlib import Path
from typing import List
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.syntax import Syntax
from tinydb import TinyDB

from .utils import move_bam
from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log

console = Console()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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


def index_bam(sample_id: str, samtools: str, threads: int, bam_output: Path) -> int:
    index_string = f"{samtools} index -@ {threads} {bam_output}"
    console.print(Panel("[bold cyan]Indexing BAM file[/bold cyan]"))
    result = run_command(index_string, f"Indexing BAM file for {sample_id}")

    if result != 0:
        console.print(f"[bold red]BAM indexing failed for sample {sample_id}[/bold red]")
        log_to_api(f"BAM indexing failed for sample {sample_id}", "ERROR", "bwa_align", sample_id,
                   str(bam_output.parent.parent))
        return 1
    return 0


def run_command(command: str, description: str) -> int:
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
        return 1

    return 0


def run_bwa(
        sample_id: str,
        fastq_files: List[Path],
        datadir: Path,
        reference: Path,
        bwa: str,
        samtools: str,
        db: TinyDB,
        threads: int = 16,
        sort_memory: str = "4G"
) -> int:
    @timer_with_db_log(db)
    def _run_bwa():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        bam_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = bam_dir / "tmp"
        temp_dir.mkdir(exist_ok=True)

        bam_output = bam_dir / f"{sample_id}.bam"
        metrics_file = bam_dir / f"{sample_id}.markdup_metrics.txt"

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
            (f"{bwa} mem -t {threads} -Y -R '{bam_header}' {reference} "
             f"{' '.join([str(file) for file in fastq_files])} > {temp_dir}/{sample_id}.sam",
             f"BWA alignment for {sample_id}"),

            (f"{samtools} view -@ {threads // 2} -bh {temp_dir}/{sample_id}.sam > {temp_dir}/{sample_id}.bam",
             f"Convert SAM to BAM for {sample_id}"),

            (f"{samtools} sort -@ {threads // 2} -n -o {temp_dir}/{sample_id}.namesort.bam {temp_dir}/{sample_id}.bam",
             f"Sort BAM by name for {sample_id}"),

            (
            f"{samtools} fixmate -@ {threads // 2} -m {temp_dir}/{sample_id}.namesort.bam {temp_dir}/{sample_id}.fixmate.bam",
            f"Fixmate for {sample_id}"),

            (f"{samtools} sort -@ {threads // 2} -m {sort_memory} -T {temp_dir}/sort_{sample_id} "
             f"-o {temp_dir}/{sample_id}.sorted.bam {temp_dir}/{sample_id}.fixmate.bam",
             f"Sort BAM by position for {sample_id}"),

            (f"{samtools} markdup -@ {threads // 2} -s {temp_dir}/{sample_id}.sorted.bam {bam_output} "
             f"2> {metrics_file}",
             f"Mark duplicates for {sample_id}"),

            (f"{samtools} index -@ {threads} {bam_output}",
             f"Index BAM for {sample_id}")
        ]

        for command, description in steps:
            if run_command(command, description) != 0:
                return 1

            output_file = command.split(">")[-1].strip().split()[0] if ">" in command else command.split()[-1]
            if output_file.endswith('.bam') or output_file.endswith('.sam'):
                file_size = Path(output_file).stat().st_size
                console.print(f"[green]File size of {output_file}: {file_size} bytes[/green]")

        shutil.rmtree(temp_dir)

        final_size = bam_output.stat().st_size
        console.print(f"[bold green]Final BAM file size: {final_size} bytes[/bold green]")

        if final_size < 1_000_000:  # Less than 1MB
            console.print("[bold red]Warning: Final BAM file is suspiciously small![/bold red]")
            return 1

        return 0

    return _run_bwa()


def print_summary(timings: dict):
    table = Table(title="BWA Alignment Process Summary")
    table.add_column("Step", style="cyan", no_wrap=True)
    table.add_column("Time (seconds)", style="magenta")

    for step, time in timings.items():
        table.add_row(step, f"{time:.2f}")

    console.print(table)


if __name__ == "__main__":
    pass