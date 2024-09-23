import subprocess
from pathlib import Path
from typing import Dict
from datetime import datetime
import psutil
import shlex
from typing import Union

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax

from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log

console = Console()

def get_gatk_version(gatk: str) -> str:
    """Get GATK version."""
    result = subprocess.run([gatk, "--version"], capture_output=True, text=True)
    return result.stdout.strip()

def haplotype_caller(datadir: Path, sample_id: str, reference: Union[str, Path], bed_file: Path, gatk: str, db: Dict, threads: int = 4, max_retries: int = 3) -> str:
    reference = Path(reference) if isinstance(reference, str) else reference

    @timer_with_db_log(db)
    def _haplotype_caller():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_vcf = vcf_dir / f"{sample_id}_GATK.vcf.gz"
        input_bam = bam_dir / f"{sample_id}.bam"
        dbsnp_file = Path("/apps/data/src/bundle/00-All.vcf.gz")

        gatk_version = get_gatk_version(gatk)
        log_to_db(db, f"Starting HaplotypeCaller for sample {sample_id} with GATK version {gatk_version}", "INFO", "GATK", sample_id, datadir.name)

        if output_vcf.exists():
            console.print(Panel(f"[yellow]VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("GATK VCF file exists", "INFO", "GATK", sample_id, datadir.name)
            log_to_db(db, f"VCF file already exists for {sample_id}", "INFO", "GATK", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting variant calling with GATK for {sample_id}[/bold blue]"))
        log_to_api("Start variant calling with GATK", "INFO", "GATK", sample_id, datadir.name)
        log_to_db(db, f"Starting variant calling with GATK for {sample_id}", "INFO", "GATK", sample_id, datadir.name)

        # Update reference dictionary
        dict_file = reference.with_suffix('.dict')
        if dict_file.exists():
            console.print(Panel(f"[green]Reference sequence dictionary already exists: {dict_file}[/green]"))
            log_to_db(db, f"Using existing reference dictionary: {dict_file}", "INFO", "GATK", sample_id, datadir.name)
        else:
            error_msg = f"Reference dictionary file not found: {dict_file}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_db(db, error_msg, "ERROR", "GATK", sample_id, datadir.name)
            log_to_api(error_msg, "ERROR", "GATK", sample_id, datadir.name)
            return "error"


        gatk_command = (
            f"{gatk} HaplotypeCaller "
            f"--input {input_bam} "
            f"--output {output_vcf} "
            f"--reference {reference} "
            f"--intervals {bed_file} "
            f"--dbsnp {dbsnp_file} "
            f"--native-pair-hmm-threads 4 "
            f"--annotation-group StandardAnnotation "
            f"--annotation StrandBiasBySample "
            f"--standard-min-confidence-threshold-for-calling 30 "
            f"--emit-ref-confidence GVCF "
            f"--create-output-variant-index true "
            f"--pcr-indel-model AGGRESSIVE "
            f"--max-alternate-alleles 3 "
            f"--contamination-fraction-to-filter 0.0 "
            f"--dont-use-soft-clipped-bases "
            f"--pairHMM LOGLESS_CACHING "
            f"--activity-profile-out {vcf_dir}/{sample_id}_activity.igv.gz "
            f"--assembly-region-out {vcf_dir}/{sample_id}_assembly.igv.gz"
        )

        console.print(Syntax(gatk_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"GATK command: {gatk_command}", "INFO", "GATK", sample_id, datadir.name)


        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"GATK HaplotypeCaller started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "GATK", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running GATK HaplotypeCaller for {sample_id}...", total=None)

                process = subprocess.Popen(shlex.split(gatk_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

                for line in process.stdout:
                    progress.update(task, advance=1)
                    console.print(f"[dim]{line.strip()}[/dim]")
                    log_to_db(db, line.strip(), "DEBUG", "GATK", sample_id, datadir.name)

                process.wait()

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"GATK HaplotypeCaller finished at {end_time}. Duration: {duration}", "INFO", "GATK", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"GATK HaplotypeCaller failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "GATK", sample_id, datadir.name)
            else:
                error_msg = f"GATK HaplotypeCaller failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "GATK", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "GATK", sample_id, datadir.name)
                return "error"

        console.print(Panel(f"[bold green]GATK variants determined for {sample_id}[/bold green]"))
        log_to_api("GATK variants determined", "INFO", "GATK", sample_id, datadir.name)
        log_to_db(db, f"GATK variants determined successfully for {sample_id}", "INFO", "GATK", sample_id, datadir.name)

        # Log file sizes and resource usage
        log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "GATK", sample_id, datadir.name)
        log_to_db(db, f"Output VCF index size: {(output_vcf.with_suffix('.vcf.gz.tbi')).stat().st_size} bytes", "INFO", "GATK", sample_id, datadir.name)
        log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "GATK", sample_id, datadir.name)

        return "success"

    return _haplotype_caller()

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass