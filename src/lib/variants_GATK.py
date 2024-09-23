import subprocess
from pathlib import Path
from typing import Dict
from datetime import datetime

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax

from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log

console = Console()

def haplotype_caller(datadir: Path, sample_id: str, reference: Path, bed_file: Path, gatk: str, db: Dict) -> str:
    """
    Function that calls GATK's HaplotypeCaller to generate a list of raw variants.

    :param datadir: The directory where the data is located.
    :param sample_id: ID of the patient/sample being analysed using GATK.
    :param reference: Reference file used in the original alignment.
    :param bed_file: BED file with regions to be analysed.
    :param gatk: GATK executable location.
    :param db: Database for logging.

    :return: Returns 'success' if the operation is successful, 'exists' if the VCF file already exists.
    """
    @timer_with_db_log(db)
    def _haplotype_caller():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_vcf = vcf_dir / f"{sample_id}_GATK.vcf.gz"
        input_bam = bam_dir / f"{sample_id}.bam"

        log_to_db(db, f"Starting HaplotypeCaller for sample {sample_id}", "INFO", "GATK", sample_id, datadir.name)

        if output_vcf.exists():
            console.print(Panel(f"[yellow]VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("GATK VCF file exists", "INFO", "GATK", sample_id, datadir.name)
            log_to_db(db, f"VCF file already exists for {sample_id}", "INFO", "GATK", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting variant calling with GATK for {sample_id}[/bold blue]"))
        log_to_api("Start variant calling with GATK", "INFO", "GATK", sample_id, datadir.name)
        log_to_db(db, f"Starting variant calling with GATK for {sample_id}", "INFO", "GATK", sample_id, datadir.name)

        gatk_command = (
            f"{gatk} HaplotypeCaller "
            f"-R {reference} "
            f"-I {input_bam} "
            f"-O {output_vcf} "
            f"-L {bed_file} "
            f"--native-pair-hmm-threads 4 "
            f"--smith-waterman FASTEST_AVAILABLE "
            f"--germline-resource /apps/data/src/bundle/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz "  
            f"--dbsnp /apps/data/src/bundle/00-All.vcf.gz "
            f"--genotyping-mode DISCOVERY "
            f"--annotation-group StandardAnnotation "
            f"--annotation-group AS_StandardAnnotation "
            f"--annotation StrandBiasBySample "
            f"-ERC GVCF "
            f"--create-output-variant-index true"
        )

        console.print(Syntax(gatk_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"GATK command: {gatk_command}", "INFO", "GATK", sample_id, datadir.name)

        start_time = datetime.now()
        log_to_db(db, f"GATK HaplotypeCaller started at {start_time}", "INFO", "GATK", sample_id, datadir.name)

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

            process = subprocess.Popen(gatk_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

            for line in process.stdout:
                progress.update(task, advance=1)
                console.print(f"[dim]{line.strip()}[/dim]")
                log_to_db(db, line.strip(), "DEBUG", "GATK", sample_id, datadir.name)

            process.wait()

        end_time = datetime.now()
        duration = end_time - start_time
        log_to_db(db, f"GATK HaplotypeCaller finished at {end_time}. Duration: {duration}", "INFO", "GATK", sample_id, datadir.name)

        if process.returncode != 0:
            error_msg = f"GATK HaplotypeCaller failed for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "GATK", sample_id, datadir.name)
            log_to_db(db, f"GATK HaplotypeCaller failed for {sample_id}. Return code: {process.returncode}", "ERROR", "GATK", sample_id, datadir.name)
            return "error"

        console.print(Panel(f"[bold green]GATK variants determined for {sample_id}[/bold green]"))
        log_to_api("GATK variants determined", "INFO", "GATK", sample_id, datadir.name)
        log_to_db(db, f"GATK variants determined successfully for {sample_id}", "INFO", "GATK", sample_id, datadir.name)

        # Log file sizes
        log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "GATK", sample_id, datadir.name)
        log_to_db(db, f"Output VCF index size: {(output_vcf.with_suffix('.vcf.gz.tbi')).stat().st_size} bytes", "INFO", "GATK", sample_id, datadir.name)

        return "success"

    return _haplotype_caller()

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass