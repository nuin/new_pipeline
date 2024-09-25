"""
.. module:: GATK_vcf
    :platform: any
    :synopsis: Module that calls GATK to compare VCF files and performs post-analysis of these variants
.. moduleauthor:: Paulo Nuin, April 2016, Updated September 2024
"""

import subprocess
from pathlib import Path
from typing import Dict
from datetime import datetime
import psutil
import shlex

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax

from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log

console = Console()

def get_gatk_version(gatk: str) -> str:
    """Get GATK version."""
    result = subprocess.run(["java", "-jar", gatk, "--version"], capture_output=True, text=True)
    return result.stdout.strip()

def vcf_comparison(datadir: Path, sample_id: str, reference: Path, gatk: str, db: Dict, max_retries: int = 3) -> str:
    """
    Function that merges the available VCFs in the sample VCF datadir

    :param datadir: Location of the BAM files
    :param sample_id: ID of the patient/sample being analysed using GATK
    :param reference: Reference file used in the original alignment
    :param gatk: GATK jar file location
    :param db: Database for logging
    :param max_retries: Maximum number of retries for GATK execution

    :return: returns 'success' if the VCF files are successfully merged, 'exists' if the merged file already exists.
    """
    @timer_with_db_log(db)
    def _vcf_comparison():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        output_vcf = vcf_dir / f"{sample_id}_merged.vcf"

        gatk_version = get_gatk_version(gatk)
        log_to_db(db, f"Starting VCF comparison for sample {sample_id} with GATK version {gatk_version}", "INFO", "GATK_vcf", sample_id, datadir.name)

        if output_vcf.exists():
            console.print(Panel(f"[yellow]Merged VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("Merged VCF file exists", "INFO", "GATK_vcf", sample_id, datadir.name)
            log_to_db(db, f"Merged VCF file already exists for {sample_id}", "INFO", "GATK_vcf", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting merge of GATK, GATK3, Freebayes, and Octopus VCFs for {sample_id}[/bold blue]"))
        log_to_api("Starting merge of GATK, GATK3, Freebayes, and Octopus VCFs", "INFO", "GATK_vcf", sample_id, datadir.name)
        log_to_db(db, f"Starting merge of GATK, GATK3, Freebayes, and Octopus VCFs for {sample_id}", "INFO", "GATK_vcf", sample_id, datadir.name)

        gatk_command = (
            f"{gatk} -T CombineVariants "
            f"-R {reference} "
            f"--variant:freebayes {vcf_dir}/{sample_id}_freebayes.vcf "
            f"--variant:gatk {vcf_dir}/{sample_id}_GATK.vcf "
            f"--variant:gatk3 {vcf_dir}/{sample_id}_GATK3.vcf "
            f"--variant:octopus {vcf_dir}/{sample_id}_octopus.vcf "
            f"-o {output_vcf} "
            f"--genotypemergeoption UNSORTED "
            f"--mergeInfoWithMaxAC "
            f"--minimumN 2"
        )

        console.print(Syntax(gatk_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"GATK command: {gatk_command}", "INFO", "GATK_vcf", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"GATK CombineVariants started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "GATK_vcf", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running GATK CombineVariants for {sample_id}...", total=None)

                process = subprocess.Popen(shlex.split(gatk_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

                for line in process.stdout:
                    progress.update(task, advance=1)
                    console.print(f"[dim]{line.strip()}[/dim]")
                    log_to_db(db, line.strip(), "DEBUG", "GATK_vcf", sample_id, datadir.name)

                process.wait()

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"GATK CombineVariants finished at {end_time}. Duration: {duration}", "INFO", "GATK_vcf", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"GATK CombineVariants failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "GATK_vcf", sample_id, datadir.name)
            else:
                error_msg = f"GATK CombineVariants failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "GATK_vcf", sample_id, datadir.name)
                return "error"

        if output_vcf.exists():
            console.print(Panel(f"[bold green]VCFs merged successfully for {sample_id}[/bold green]"))
            log_to_api("VCFs merged successfully", "INFO", "GATK_vcf", sample_id, datadir.name)
            log_to_db(db, f"VCFs merged successfully for {sample_id}", "INFO", "GATK_vcf", sample_id, datadir.name)

            # Log file size and resource usage
            log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "GATK_vcf", sample_id, datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "GATK_vcf", sample_id, datadir.name)

            return "success"
        else:
            error_msg = f"GATK CombineVariants completed but output VCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
            return "error"

    return _vcf_comparison()

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass