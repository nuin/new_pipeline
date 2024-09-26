"""
.. module:: vep_ann
    :platform: any
    :synopsis: This module calls Ensembl Variant Effect Predictor (VEP) to annotate variants from the merged VCF file
.. moduleauthor:: Assistant, September 2024
"""

import subprocess
from pathlib import Path
from typing import Dict, Optional
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

def get_vep_version(vep: str) -> str:
    """Get VEP version."""
    result = subprocess.run([vep, "--help"], capture_output=True, text=True)
    for line in result.stdout.split('\n'):
        if "ensembl-vep" in line:
            return line.strip()
    return "Unknown"

def annotate_merged(
    sample_id: str,
    datadir: Path,
    vep: str,
    cache_dir: Path,
    fasta: Path,
    db: Dict,
    max_retries: int = 3
) -> str:
    @timer_with_db_log(db)
    def _annotate_merged():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        input_vcf = vcf_dir / f"{sample_id}_merged.vcf"
        output_vcf = vcf_dir / f"{sample_id}_merged.vep.vcf"

        vep_version = get_vep_version(vep)
        log_to_db(db, f"Starting VEP annotation for sample {sample_id} with VEP version {vep_version}", "INFO", "vep_ann", sample_id, datadir.name)

        if output_vcf.exists():
            console.print(Panel(f"[yellow]VEP annotated VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("VEP annotated VCF file exists", "INFO", "vep_ann", sample_id, datadir.name)
            log_to_db(db, f"VEP annotated VCF file already exists for {sample_id}", "INFO", "vep_ann", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting VEP annotation for {sample_id}[/bold blue]"))
        log_to_api("Start VEP annotation", "INFO", "vep_ann", sample_id, datadir.name)
        log_to_db(db, f"Starting VEP annotation for {sample_id}", "INFO", "vep_ann", sample_id, datadir.name)

        vep_cmd = [
            vep,
            "-i", str(input_vcf),
            "-o", str(output_vcf),
            "--offline",
            "--cache",
            "--dir_cache", str(cache_dir),
            "--assembly", "GRCh37",
            "--species", "homo_sapiens",
            "--fasta", str(fasta),
            "--format", "vcf",
            "--vcf",
            "--force_overwrite",
            "--hgvs"
        ]

        console.print(Syntax(" ".join(map(str, vep_cmd)), "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"VEP command: {' '.join(map(str, vep_cmd))}", "INFO", "vep_ann", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"VEP annotation started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "vep_ann", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running VEP annotation for {sample_id}...", total=None)

                process = subprocess.Popen(vep_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)

                stdout, stderr = process.communicate()

                for line in stdout.splitlines():
                    progress.update(task, advance=1)
                    console.print(f"[dim]{line.strip()}[/dim]")
                    log_to_db(db, line.strip(), "DEBUG", "vep_ann", sample_id, datadir.name)

                if stderr:
                    for line in stderr.splitlines():
                        console.print(f"[yellow]{line.strip()}[/yellow]")
                        log_to_db(db, line.strip(), "WARNING", "vep_ann", sample_id, datadir.name)

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"VEP annotation finished at {end_time}. Duration: {duration}", "INFO", "vep_ann", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"VEP annotation failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "vep_ann", sample_id, datadir.name)
            else:
                error_msg = f"VEP annotation failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "vep_ann", sample_id, datadir.name)
                return "error"

        if output_vcf.exists():
            console.print(Panel(f"[bold green]VEP annotation completed for {sample_id}[/bold green]"))
            log_to_api("VEP annotation completed", "INFO", "vep_ann", sample_id, datadir.name)
            log_to_db(db, f"VEP annotation completed successfully for {sample_id}", "INFO", "vep_ann", sample_id, datadir.name)

            # Log file size and resource usage
            log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "vep_ann", sample_id, datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "vep_ann", sample_id, datadir.name)

            return "success"
        else:
            error_msg = f"VEP annotation completed but output VCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
            return "error"

    return _annotate_merged()

if __name__ == "__main__":
    # Example usage (for testing purposes)
    from pathlib import Path
    from tinydb import TinyDB

    sample_id = "test_sample"
    datadir = Path("/path/to/data")
    vep = "/path/to/vep"
    cache_dir = Path("/path/to/vep/cache")
    fasta = Path("/path/to/reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz")
    db = TinyDB("/path/to/database.json")

    result = annotate_merged(sample_id, datadir, vep, cache_dir, fasta, db)
    print(f"Annotation result: {result}")