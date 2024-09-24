"""
.. module:: variants_freebayes
    :platform: Any
    :synopsis: Module that generates variants by calling Freebayes, sorts and edits the resulting VCF
.. moduleauthor:: Paulo Nuin, January 2016, Updated September 2024
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

def get_freebayes_version(freebayes: str) -> str:
    """Get Freebayes version."""
    result = subprocess.run([freebayes, "--version"], capture_output=True, text=True)
    return result.stdout.strip()

def freebayes_caller(datadir: Path, sample_id: str, reference: Path, bed_file: Path, freebayes: str, db: Dict, max_retries: int = 3) -> str:
    @timer_with_db_log(db)
    def _freebayes_caller():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_vcf = vcf_dir / f"{sample_id}_freebayes.vcf"
        input_bam = bam_dir / f"{sample_id}.bam"

        freebayes_version = get_freebayes_version(freebayes)
        log_to_db(db, f"Starting Freebayes for sample {sample_id} with Freebayes version {freebayes_version}", "INFO", "Freebayes", sample_id, datadir.name)

        if output_vcf.exists():
            console.print(Panel(f"[yellow]Freebayes VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("Freebayes VCF file exists", "INFO", "Freebayes", sample_id, datadir.name)
            log_to_db(db, f"Freebayes VCF file already exists for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting variant calling with Freebayes for {sample_id}[/bold blue]"))
        log_to_api("Start variant calling with Freebayes", "INFO", "Freebayes", sample_id, datadir.name)
        log_to_db(db, f"Starting variant calling with Freebayes for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)

        freebayes_command = (
            f"{freebayes} "
            f"-f {reference} "  # Reference genome
            f"-b {input_bam} "  # Input BAM file
            f"-t {bed_file} "   # Target regions
            f"-v {output_vcf} " # Output VCF file
            f"--min-alternate-fraction 0.2 "  # Minimum alternate allele fraction
            f"--min-alternate-count 2 "       # Minimum number of alternate allele observations
            f"--min-base-quality 20 "         # Minimum base quality to be considered
            f"--min-mapping-quality 20 "      # Minimum mapping quality to be considered
            f"--min-coverage 10 "             # Minimum coverage to consider a site
            f"--genotype-qualities "          # Calculate genotype qualities (GQ field)
            f"--strict-vcf "                  # Use strict VCF format (for compatibility)
            f"--pooled-continuous "           # Output all alleles which pass input filters
            f"--use-best-n-alleles 4 "        # Use only best N alleles, ranked by sum of supporting quality scores
            f"--haplotype-length 50 "         # Allow haplotype calls with contiguous embedded matches of up to this length
            f"--report-genotype-likelihood-max "  # Report genotypes using the maximum-likelihood estimate
            f"--theta 0.001 "                 # The expected mutation rate or pairwise nucleotide diversity among the population
            f"--ploidy 2 "                    # Specify the ploidy of the sample
        )

        console.print(Syntax(freebayes_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"Freebayes command: {freebayes_command}", "INFO", "Freebayes", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"Freebayes started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "Freebayes", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running Freebayes for {sample_id}...", total=None)

                process = subprocess.Popen(shlex.split(freebayes_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

                for line in process.stdout:
                    progress.update(task, advance=1)
                    console.print(f"[dim]{line.strip()}[/dim]")
                    log_to_db(db, line.strip(), "DEBUG", "Freebayes", sample_id, datadir.name)

                process.wait()

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"Freebayes finished at {end_time}. Duration: {duration}", "INFO", "Freebayes", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"Freebayes failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "Freebayes", sample_id, datadir.name)
            else:
                error_msg = f"Freebayes failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "Freebayes", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "Freebayes", sample_id, datadir.name)
                return "error"

        if output_vcf.exists() and output_vcf.stat().st_size > 0:
            console.print(Panel(f"[bold green]Freebayes variants determined for {sample_id}[/bold green]"))
            log_to_api("Freebayes variants determined", "INFO", "Freebayes", sample_id, datadir.name)
            log_to_db(db, f"Freebayes variants determined successfully for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)

            # Log file size and resource usage
            log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "Freebayes", sample_id, datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "Freebayes", sample_id, datadir.name)

            return "success"
        else:
            error_msg = f"Freebayes completed but output VCF not found or empty for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "Freebayes", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "Freebayes", sample_id, datadir.name)
            return "error"

    return _freebayes_caller()

def process_freebayes_vcf(datadir: Path, sample_id: str, reference: Path, db: Dict) -> str:
    @timer_with_db_log(db)
    def _process_freebayes_vcf():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        input_vcf = vcf_dir / f"{sample_id}_freebayes.vcf"
        sorted_vcf = vcf_dir / f"{sample_id}_freebayes.sorted.vcf"
        reference_dict = reference.with_suffix('.dict')

        if sorted_vcf.exists():
            console.print(Panel(f"[yellow]Sorted Freebayes VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("Sorted Freebayes VCF file exists", "INFO", "Freebayes", sample_id, datadir.name)
            log_to_db(db, f"Sorted Freebayes VCF file already exists for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)
            return "exists"

        # Step 1: Sort VCF
        try:
            console.print(Panel(f"[bold blue]Sorting Freebayes VCF for {sample_id}[/bold blue]"))
            log_to_api(f"Start sorting Freebayes VCF for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)
            log_to_db(db, f"Start sorting Freebayes VCF for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)

            # Read the reference dictionary to get the correct chromosome order
            chrom_order = []
            with reference_dict.open() as f:
                for line in f:
                    if line.startswith("@SQ"):
                        chrom = line.split("\t")[1].split(":")[1]
                        chrom_order.append(chrom)

            # Read the VCF file
            headers = []
            variants = []
            with input_vcf.open() as f:
                for line in f:
                    if line.startswith("#"):
                        if not line.startswith("##contig=<ID"):  # Skip contig lines
                            headers.append(line)
                    else:
                        variants.append(line.split("\t"))

            # Sort the variants
            sorted_variants = sorted(variants, key=lambda x: (chrom_order.index(x[0]) if x[0] in chrom_order else len(chrom_order), int(x[1])))

            # Write the sorted VCF
            with sorted_vcf.open("w") as f:
                f.writelines(headers)
                f.writelines("\t".join(variant) for variant in sorted_variants)

            console.print(Panel(f"[bold green]Freebayes VCF sorted successfully for {sample_id}[/bold green]"))
            log_to_api(f"Freebayes VCF sorted successfully for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)
            log_to_db(db, f"Freebayes VCF sorted successfully for {sample_id}", "INFO", "Freebayes", sample_id, datadir.name)

            # Log file sizes
            log_to_db(db, f"Input VCF size: {input_vcf.stat().st_size} bytes", "INFO", "Freebayes", sample_id, datadir.name)
            log_to_db(db, f"Sorted VCF size: {sorted_vcf.stat().st_size} bytes", "INFO", "Freebayes", sample_id, datadir.name)

            return "success"
        except Exception as e:
            error_msg = f"Error processing Freebayes VCF for {sample_id}: {str(e)}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "Freebayes", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "Freebayes", sample_id, datadir.name)
            return "error"

    return _process_freebayes_vcf()



if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass