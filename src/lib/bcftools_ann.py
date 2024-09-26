"""
.. module:: bcftools_ann
    :platform: any
    :synopsis: This module calls bcftools csq to annotate variants from the merged VCF file
.. moduleauthor:: Assistant, September 2024
"""

import subprocess
from pathlib import Path
import tempfile
import os
from typing import Dict, Optional
from datetime import datetime

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax

from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log


import psutil
import shlex



console = Console()

def preprocess_vcf(input_vcf: Path, output_vcf: Path, bcftools: str, db: Dict, sample_id: str) -> bool:
    cmd = [bcftools, "view", "-O", "v", "-o", str(output_vcf), str(input_vcf)]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        log_to_db(db, f"Preprocessed VCF: {output_vcf}", "INFO", "bcftools_ann", sample_id, input_vcf.parent.parent.name)
        if result.stderr:
            log_to_db(db, f"Preprocessing warnings: {result.stderr}", "WARNING", "bcftools_ann", sample_id, input_vcf.parent.parent.name)
        return True
    except subprocess.CalledProcessError as e:
        console.print(Panel(f"[bold red]Error preprocessing VCF: {e.stderr}[/bold red]"))
        log_to_db(db, f"Error preprocessing VCF: {e.stderr}", "ERROR", "bcftools_ann", sample_id, input_vcf.parent.parent.name)
        return False
def get_bcftools_version(bcftools: str) -> str:
    """Get bcftools version."""
    result = subprocess.run([bcftools, "--version"], capture_output=True, text=True)
    return result.stdout.split('\n')[0]

def annotate_merged(
    sample_id: str,
    datadir: Path,
    bcftools: str,
    reference: Path,
    gff: Path,
    db: Dict,
    transcript_list: Optional[Path] = None,
    max_retries: int = 3
) -> str:
    @timer_with_db_log(db)
    def _annotate_merged():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        input_vcf = vcf_dir / f"{sample_id}_merged.vcf"
        output_vcf = vcf_dir / f"{sample_id}_merged.csq.vcf"
        final_vcf = vcf_dir / f"{sample_id}_merged.csq.filtered.vcf"

        bcftools_version = get_bcftools_version(bcftools)
        log_to_db(db, f"Starting bcftools csq for sample {sample_id} with bcftools version {bcftools_version}", "INFO", "bcftools_ann", sample_id, datadir.name)

        if final_vcf.exists():
            console.print(Panel(f"[yellow]Annotated VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("Annotated VCF file exists", "INFO", "bcftools_ann", sample_id, datadir.name)
            log_to_db(db, f"Annotated VCF file already exists for {sample_id}", "INFO", "bcftools_ann", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting bcftools csq annotation for {sample_id}[/bold blue]"))
        log_to_api("Start bcftools csq annotation", "INFO", "bcftools_ann", sample_id, datadir.name)
        log_to_db(db, f"Starting bcftools csq annotation for {sample_id}", "INFO", "bcftools_ann", sample_id, datadir.name)

        bcftools_cmd = [
            bcftools, "csq",
            "-f", str(reference),
            "-g", str(gff),
            "--local-csq",
            "--ncsq", "20",
            "-Ov",
            "-o", str(output_vcf),
            str(input_vcf)
        ]

        console.print(Syntax(" ".join(map(str, bcftools_cmd)), "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"bcftools command: {' '.join(map(str, bcftools_cmd))}", "INFO", "bcftools_ann", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"bcftools csq started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "bcftools_ann", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running bcftools csq for {sample_id}...", total=None)

                process = subprocess.Popen(bcftools_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)

                stdout, stderr = process.communicate()

                for line in stdout.splitlines():
                    progress.update(task, advance=1)
                    console.print(f"[dim]{line.strip()}[/dim]")
                    log_to_db(db, line.strip(), "DEBUG", "bcftools_ann", sample_id, datadir.name)

                if stderr:
                    for line in stderr.splitlines():
                        console.print(f"[yellow]{line.strip()}[/yellow]")
                        log_to_db(db, line.strip(), "WARNING", "bcftools_ann", sample_id, datadir.name)

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"bcftools csq finished at {end_time}. Duration: {duration}", "INFO", "bcftools_ann", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"bcftools csq failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "bcftools_ann", sample_id, datadir.name)
            else:
                error_msg = f"bcftools csq failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "bcftools_ann", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "bcftools_ann", sample_id, datadir.name)
                return "error"

        if output_vcf.exists():
            if transcript_list:
                console.print(Panel(f"[bold blue]Filtering VCF for specific transcripts for {sample_id}[/bold blue]"))
                log_to_api("Filtering VCF for specific transcripts", "INFO", "bcftools_ann", sample_id, datadir.name)
                log_to_db(db, f"Filtering VCF for specific transcripts for {sample_id}", "INFO", "bcftools_ann", sample_id, datadir.name)

                try:
                    # Read the list of desired transcripts
                    with open(transcript_list, 'r') as f:
                        desired_transcripts = set(line.strip() for line in f)

                    # Filter the CSQ output
                    with open(output_vcf, 'r') as infile, open(final_vcf, 'w') as outfile:
                        for line in infile:
                            if line.startswith('#'):
                                outfile.write(line)
                                continue
                            
                            fields = line.strip().split('\t')
                            info = fields[7]
                            csq_field = next((f for f in info.split(';') if f.startswith('BCSQ=')), None)
                            
                            if csq_field:
                                csq_parts = csq_field.split('|')
                                transcript = csq_parts[4] if len(csq_parts) > 4 else None
                                
                                if transcript in desired_transcripts:
                                    outfile.write(line)
                            else:
                                outfile.write(line)

                    console.print(Panel(f"[bold green]VCF filtered for specific transcripts for {sample_id}[/bold green]"))
                    log_to_api("VCF filtered for specific transcripts", "INFO", "bcftools_ann", sample_id, datadir.name)
                    log_to_db(db, f"VCF filtered for specific transcripts for {sample_id}", "INFO", "bcftools_ann", sample_id, datadir.name)

                    # Replace the original CSQ output with the filtered one
                    final_vcf.replace(output_vcf)
                except Exception as e:
                    error_msg = f"Error filtering VCF for specific transcripts: {str(e)}"
                    console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                    log_to_api(error_msg, "ERROR", "bcftools_ann", sample_id, datadir.name)
                    log_to_db(db, error_msg, "ERROR", "bcftools_ann", sample_id, datadir.name)
                    return "error"
            else:
                final_vcf = output_vcf

            console.print(Panel(f"[bold green]bcftools csq annotation completed for {sample_id}[/bold green]"))
            log_to_api("bcftools csq annotation completed", "INFO", "bcftools_ann", sample_id, datadir.name)
            log_to_db(db, f"bcftools csq annotation completed successfully for {sample_id}", "INFO", "bcftools_ann", sample_id, datadir.name)

            # Log file size and resource usage
            log_to_db(db, f"Final VCF size: {final_vcf.stat().st_size} bytes", "INFO", "bcftools_ann", sample_id, datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "bcftools_ann", sample_id, datadir.name)

            return "success"
        else:
            error_msg = f"bcftools csq completed but output VCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "bcftools_ann", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "bcftools_ann", sample_id, datadir.name)
            return "error"

    return _annotate_merged()



if __name__ == "__main__":
    # Example usage (for testing purposes)
    from pathlib import Path
    from tinydb import TinyDB

    sample_id = "test_sample"
    datadir = Path("/path/to/data")
    bcftools = "/path/to/bcftools"
    reference = Path("/path/to/reference.fa")
    gff = Path("/path/to/annotation.gff3.gz")
    db = TinyDB("/path/to/database.json")

    result = annotate_merged(sample_id, datadir, bcftools, reference, gff, db)
    print(f"Annotation result: {result}")