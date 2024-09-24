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

def get_octopus_version(octopus: str) -> str:
    """Get Octopus version."""
    try:
        result = subprocess.run([octopus, "--version"], capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return "Unknown"

def octopus_caller(datadir: Path, sample_id: str, reference: Path, bed_file: Path, octopus: str, db: Dict, threads: int = 16, max_retries: int = 3) -> str:
    @timer_with_db_log(db)
    def _octopus_caller():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_vcf = vcf_dir / f"{sample_id}_octopus.vcf"
        input_bam = bam_dir / f"{sample_id}.bam"

        octopus_version = get_octopus_version(octopus)
        log_to_db(db, f"Starting Octopus for sample {sample_id} with Octopus version {octopus_version}", "INFO", "Octopus", sample_id, datadir.name)

        if output_vcf.exists():
            console.print(Panel(f"[yellow]Octopus VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("Octopus VCF file exists", "INFO", "Octopus", sample_id, datadir.name)
            log_to_db(db, f"Octopus VCF file already exists for {sample_id}", "INFO", "Octopus", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting germline variant calling with Octopus for {sample_id}[/bold blue]"))
        log_to_api("Start germline variant calling with Octopus", "INFO", "Octopus", sample_id, datadir.name)
        log_to_db(db, f"Starting germline variant calling with Octopus for {sample_id}", "INFO", "Octopus", sample_id, datadir.name)

        octopus_command = (
            f"{octopus} "
            f"-R {reference} "
            f"-I {input_bam} "
            f"--regions-file {bed_file} "
            f"--threads {threads} "
            f"-o {output_vcf} "
            f"--min-variant-posterior 0.01 "
            f"--annotations AD DP ADP GQ GT MQ AF AC AN SB BQ "
            f"ABP ADP ADRP ARF AOR CYC DAD DC DCP DMP DPF ED FEAD FRF GC "
            f"MP MRC PP QD QUAL REB RSB RTB SD SF SHC SMQ"
        )

        console.print(Syntax(octopus_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"Octopus command: {octopus_command}", "INFO", "Octopus", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"Octopus started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "Octopus", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running Octopus for {sample_id}...", total=None)

                process = subprocess.Popen(shlex.split(octopus_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=False)

                while True:
                    output = process.stdout.readline()
                    if output == b'' and process.poll() is not None:
                        break
                    if output:
                        try:
                            line = output.decode('utf-8').strip()
                        except UnicodeDecodeError:
                            line = output.decode('utf-8', errors='replace').strip()
                        progress.update(task, advance=1)
                        console.print(f"[dim]{line}[/dim]")
                        log_to_db(db, line, "DEBUG", "Octopus", sample_id, datadir.name)

                process.wait()

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"Octopus finished at {end_time}. Duration: {duration}", "INFO", "Octopus", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"Octopus failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "Octopus", sample_id, datadir.name)
            else:
                error_msg = f"Octopus failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "Octopus", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "Octopus", sample_id, datadir.name)
                return "error"

        if output_vcf.exists():
            console.print(Panel(f"[bold green]Octopus germline variants determined for {sample_id}[/bold green]"))
            log_to_api("Octopus germline variants determined", "INFO", "Octopus", sample_id, datadir.name)
            log_to_db(db, f"Octopus germline variants determined successfully for {sample_id}", "INFO", "Octopus", sample_id, datadir.name)

            # Log file size and resource usage
            log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "Octopus", sample_id, datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "Octopus", sample_id, datadir.name)

            # Change VCF version
            change_vcf_version(output_vcf)

            return "success"
        else:
            error_msg = f"Octopus completed but output VCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "Octopus", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "Octopus", sample_id, datadir.name)
            return "error"

    return _octopus_caller()

def change_vcf_version(vcf_file: Path) -> None:
    """
    Function that changes the VCF version from 4.3 to 4.2 so it can be used in GATK

    :param vcf_file: VCF file location
    :type vcf_file: Path

    :return: None
    """
    console.print(Panel(f"[bold blue]Changing VCF version for compatibility: {vcf_file}[/bold blue]"))

    temp_file = vcf_file.with_suffix('.vcf.temp')

    try:
        with vcf_file.open('r') as infile, temp_file.open('w') as outfile:
            for line in infile:
                if line.startswith('##fileformat='):
                    outfile.write('##fileformat=VCFv4.2\n')
                else:
                    outfile.write(line)

        # Replace the original file with the modified one
        temp_file.replace(vcf_file)

        console.print(Panel(f"[bold green]Successfully changed VCF version to 4.2 for {vcf_file}[/bold green]"))
    except Exception as e:
        console.print(Panel(f"[bold red]Error changing VCF version: {str(e)}[/bold red]"))
    finally:
        if temp_file.exists():
            temp_file.unlink()

if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass