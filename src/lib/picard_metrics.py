"""
.. module:: picard_metrics
    :platform: any
    :synopsis: Module that provides parsing capabilities for multiple types of Picard output
.. moduleauthor:: Paulo Nuin, January (modified May) 2016, Updated September 2024
"""

import subprocess
from pathlib import Path
from typing import Dict

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax

from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log

console = Console()

def run_picard_command(command: str, description: str, sample_id: str, datadir: Path, db: Dict) -> int:
    console.print(Panel(f"[bold blue]{description}[/bold blue]"))
    console.print(Syntax(command, "bash", theme="monokai", line_numbers=True))

    log_to_db(db, f"Running Picard command: {command}", "INFO", "picard_metrics", sample_id, datadir.name)

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

        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

        for line in process.stdout:
            progress.update(task, advance=1)
            console.print(f"[dim]{line.strip()}[/dim]")
            log_to_db(db, line.strip(), "DEBUG", "picard_metrics", sample_id, datadir.name)

    process.wait()

    if process.returncode != 0:
        error_msg = f"Error in {description}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_api(error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
        log_to_db(db, error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
        return 1

    return 0

def get_yield(sample_id: str, datadir: Path, picard: str, db: Dict) -> str:
    @timer_with_db_log(db)
    def _get_yield():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        output_file = metrics_dir / f"{sample_id}.yield.out"

        if output_file.exists():
            console.print(Panel(f"[yellow]Picard CollectQualityYieldMetrics file exists for {sample_id}[/yellow]"))
            log_to_api("Picard CollectQualityYieldMetrics file exists", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"Picard CollectQualityYieldMetrics file exists for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "exists"

        picard_command = (
            f"{picard} CollectQualityYieldMetrics "
            f"I={bam_dir}/{sample_id}.bam "
            f"O={output_file}"
        )

        result = run_picard_command(picard_command, f"Generating Picard CollectQualityYieldMetrics for {sample_id}", sample_id, datadir, db)

        if result == 0:
            console.print(Panel(f"[bold green]CollectQualityYieldMetrics file created for {sample_id}[/bold green]"))
            log_to_api("CollectQualityYieldMetrics file created", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"CollectQualityYieldMetrics file created for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "success"
        else:
            return "error"

    return _get_yield()


def get_hs_metrics(sample_id: str, datadir: Path, reference: Path, bait_file: Path, picard: str, db: Dict, panel: str = "full") -> str:
    @timer_with_db_log(db)
    def _get_hs_metrics():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"

        if panel == "full":
            output_file = metrics_dir / f"{sample_id}.hs_metrics.out"
        else:
            output_file = metrics_dir / f"{sample_id}.hs_metrics.panel.out"
            bait_file = bait_file.with_suffix('.picard.bed')

        if output_file.exists():
            console.print(Panel(f"[yellow]Picard CollectHsMetrics {panel} file exists for {sample_id}[/yellow]"))
            log_to_api(f"Picard CollectHsMetrics {panel} file exists", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"Picard CollectHsMetrics {panel} file exists for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "exists"

        # Ensure bait_file is a Path object and exists
        bait_file = Path(bait_file)
        if not bait_file.exists():
            error_msg = f"Bait file does not exist: {bait_file}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
            return "error"

        picard_command = (
            f"{picard} CollectHsMetrics "
            f"I={bam_dir}/{sample_id}.bam "
            f"O={output_file} "
            f"R={reference} "
            f"BAIT_INTERVALS={bait_file} "
            f"TARGET_INTERVALS={bait_file}"
        )

        result = run_picard_command(picard_command, f"Generating Picard CollectHsMetrics {panel} for {sample_id}", sample_id, datadir, db)

        if result == 0:
            console.print(Panel(f"[bold green]CollectHsMetrics {panel} file created for {sample_id}[/bold green]"))
            log_to_api(f"CollectHsMetrics {panel} file created", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"CollectHsMetrics {panel} file created for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "success"
        else:
            return "error"

    return _get_hs_metrics()


def get_align_summary(sample_id: str, datadir: Path, reference: Path, picard: str, db: Dict) -> str:
    @timer_with_db_log(db)
    def _get_align_summary():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        output_file = metrics_dir / f"{sample_id}.align_metrics.out"

        if output_file.exists():
            console.print(Panel(f"[yellow]Picard AlignmentSummaryMetrics file exists for {sample_id}[/yellow]"))
            log_to_api("Picard AlignmentSummaryMetrics file exists", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"Picard AlignmentSummaryMetrics file exists for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "exists"

        picard_command = (
            f"{picard} CollectAlignmentSummaryMetrics "
            f"I={bam_dir}/{sample_id}.bam "
            f"O={output_file} "
            f"R={reference}"
        )

        result = run_picard_command(picard_command, f"Generating Picard CollectAlignmentSummaryMetrics for {sample_id}", sample_id, datadir, db)

        if result == 0:
            console.print(Panel(f"[bold green]CollectAlignmentSummaryMetrics file created for {sample_id}[/bold green]"))
            log_to_api("CollectAlignmentSummaryMetrics file created", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"CollectAlignmentSummaryMetrics file created for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "success"
        else:
            return "error"

    return _get_align_summary()

def get_call_metrics(sample_id: str, datadir: Path, vcf_file: Path, picard: str, db: Dict) -> str:
    @timer_with_db_log(db)
    def _get_call_metrics():
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        output_file = metrics_dir / f"{sample_id}.call_metrics.out"
        input_vcf = vcf_dir / f"{sample_id}_merged.vcf"

        if (output_file.with_suffix('.variant_calling_detail_metrics')).exists():
            console.print(Panel(f"[yellow]Picard CollectVariantCallingMetrics file exists for {sample_id}[/yellow]"))
            log_to_api("Picard CollectVariantCallingMetrics file exists", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"Picard CollectVariantCallingMetrics file exists for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "exists"

        picard_command = (
            f"{picard} CollectVariantCallingMetrics "
            f"I={input_vcf} "
            f"O={output_file} "
            f"DBSNP={vcf_file} "
            f"QUIET=true"
        )

        result = run_picard_command(picard_command, f"Generating Picard CollectVariantCallingMetrics for {sample_id}", sample_id, datadir, db)

        if result == 0:
            console.print(Panel(f"[bold green]CollectVariantCallingMetrics file created for {sample_id}[/bold green]"))
            log_to_api("CollectVariantCallingMetrics file created", "INFO", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, f"CollectVariantCallingMetrics file created for {sample_id}", "INFO", "picard_metrics", sample_id, datadir.name)
            return "success"
        else:
            return "error"

    return _get_call_metrics()

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass