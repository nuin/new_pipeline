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

def get_gatk3_version(gatk: str) -> str:
    """Get GATK3 version."""
    result = subprocess.run(["java", "-jar", gatk, "--version"], capture_output=True, text=True)
    return result.stdout.strip()

def haplotype_caller(datadir: Path, sample_id: str, reference: Path, bed_file: Path, gatk: str, db: Dict, threads: int = 16, max_retries: int = 3) -> str:
    @timer_with_db_log(db)
    def _haplotype_caller():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_vcf = vcf_dir / f"{sample_id}_GATK3.vcf"
        input_bam = bam_dir / f"{sample_id}.bam"
        dbsnp_file = Path("/apps/data/src/bundle/00-All.vcf.gz")

        gatk_version = get_gatk3_version(gatk)
        log_to_db(db, f"Starting HaplotypeCaller for sample {sample_id} with GATK3 version {gatk_version}", "INFO", "GATK3", sample_id, datadir.name)

        if output_vcf.exists():
            console.print(Panel(f"[yellow]VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("GATK3 VCF file exists", "INFO", "GATK3", sample_id, datadir.name)
            log_to_db(db, f"VCF file already exists for {sample_id}", "INFO", "GATK3", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting variant calling with GATK3 for {sample_id}[/bold blue]"))
        log_to_api("Start variant calling with GATK3", "INFO", "GATK3", sample_id, datadir.name)
        log_to_db(db, f"Starting variant calling with GATK3 for {sample_id}", "INFO", "GATK3", sample_id, datadir.name)

        gatk_command = (
            f"java -jar {gatk} -T HaplotypeCaller "
            f"-R {reference} "
            f"-I {input_bam} "
            f"-o {output_vcf} "
            f"-L {bed_file} "
            f"-ip 2 "
            f"-A StrandBiasBySample "
            f"-A BaseQualityRankSumTest "
            f"-A MappingQualityRankSumTest "
            f"-A ReadPosRankSumTest "
            f"-A FisherStrand "
            f"-A QualByDepth "
            f"-A Coverage "
            f"--dbsnp {dbsnp_file} "
            f"-stand_call_conf 30 "
            f"-pairHMM VECTOR_LOGLESS_CACHING "
            f"-nct {threads} "
            f"--max_alternate_alleles 3 "
            f"--contamination_fraction_to_filter 0.0 "
            f"--pcr_indel_model AGGRESSIVE "
            f"--dontUseSoftClippedBases"
        )

        console.print(Syntax(gatk_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"GATK3 command: {gatk_command}", "INFO", "GATK3", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"GATK3 HaplotypeCaller started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "GATK3", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running GATK3 HaplotypeCaller for {sample_id}...", total=None)

                process = subprocess.Popen(shlex.split(gatk_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

                for line in process.stdout:
                    progress.update(task, advance=1)
                    console.print(f"[dim]{line.strip()}[/dim]")
                    log_to_db(db, line.strip(), "DEBUG", "GATK3", sample_id, datadir.name)

                process.wait()

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"GATK3 HaplotypeCaller finished at {end_time}. Duration: {duration}", "INFO", "GATK3", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"GATK3 HaplotypeCaller failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "GATK3", sample_id, datadir.name)
            else:
                error_msg = f"GATK3 HaplotypeCaller failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "GATK3", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "GATK3", sample_id, datadir.name)
                return "error"

        if output_vcf.exists():
            console.print(Panel(f"[bold green]GATK3 variants determined for {sample_id}[/bold green]"))
            log_to_api("GATK3 variants determined", "INFO", "GATK3", sample_id, datadir.name)
            log_to_db(db, f"GATK3 variants determined successfully for {sample_id}", "INFO", "GATK3", sample_id, datadir.name)

            # Log file size and resource usage
            log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "GATK3", sample_id, datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "GATK3", sample_id, datadir.name)

            return "success"
        else:
            error_msg = f"GATK3 HaplotypeCaller completed but output VCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "GATK3", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "GATK3", sample_id, datadir.name)
            return "error"

    return _haplotype_caller()

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass