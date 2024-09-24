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

# Replace 'docker_user' with the actual username that has Docker permissions
DOCKER_USER = 'docker_user'


def run_as_docker_user(command: str) -> subprocess.CompletedProcess:
    """Run a command as the Docker user."""
    full_command = f"su - {DOCKER_USER} -c '{command}'"
    return subprocess.run(shlex.split(full_command), capture_output=True, text=True)


def get_deepvariant_version() -> str:
    """Get DeepVariant version from Docker image."""
    command = "docker run google/deepvariant:1.6.1 /opt/deepvariant/bin/run_deepvariant --version"
    result = run_as_docker_user(command)
    return result.stdout.strip()


def deepvariant_caller(datadir: Path, sample_id: str, reference: Path, bed_file: Path, db: Dict, threads: int = 16,
                       max_retries: int = 3) -> str:
    @timer_with_db_log(db)
    def _deepvariant_caller():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_vcf = vcf_dir / f"{sample_id}_DeepVariant.vcf.gz"
        output_gvcf = vcf_dir / f"{sample_id}_DeepVariant.g.vcf.gz"
        input_bam = bam_dir / f"{sample_id}.bam"
        intermediate_results_dir = vcf_dir / "deepvariant_intermediate"

        deepvariant_version = get_deepvariant_version()
        log_to_db(db, f"Starting DeepVariant for sample {sample_id} with DeepVariant version {deepvariant_version}",
                  "INFO", "DeepVariant", sample_id, datadir.name)

        if output_vcf.exists() and output_gvcf.exists():
            console.print(Panel(f"[yellow]DeepVariant VCF and gVCF files already exist for {sample_id}[/yellow]"))
            log_to_api("DeepVariant VCF and gVCF files exist", "INFO", "DeepVariant", sample_id, datadir.name)
            log_to_db(db, f"DeepVariant VCF and gVCF files already exist for {sample_id}", "INFO", "DeepVariant",
                      sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting variant calling with DeepVariant for {sample_id}[/bold blue]"))
        log_to_api("Start variant calling with DeepVariant", "INFO", "DeepVariant", sample_id, datadir.name)
        log_to_db(db, f"Starting variant calling with DeepVariant for {sample_id}", "INFO", "DeepVariant", sample_id,
                  datadir.name)

        # Ensure output directories exist
        vcf_dir.mkdir(parents=True, exist_ok=True)
        intermediate_results_dir.mkdir(parents=True, exist_ok=True)

        deepvariant_command = (
            f"docker run "
            f"-v {datadir}:/input "
            f"-v {vcf_dir}:/output "
            f"google/deepvariant:1.6.1 "
            f"/opt/deepvariant/bin/run_deepvariant "
            f"--model_type=WGS "  # Adjust this based on your data type (WGS, WES, PACBIO, etc.)
            f"--ref=/input/{reference.relative_to(datadir)} "
            f"--reads=/input/{input_bam.relative_to(datadir)} "
            f"--regions=/input/{bed_file.relative_to(datadir)} "
            f"--output_vcf=/output/{output_vcf.name} "
            f"--output_gvcf=/output/{output_gvcf.name} "
            f"--intermediate_results_dir=/output/{intermediate_results_dir.name} "
            f"--num_shards={threads}"
        )

        console.print(Syntax(deepvariant_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"DeepVariant command: {deepvariant_command}", "INFO", "DeepVariant", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"DeepVariant started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO",
                      "DeepVariant", sample_id, datadir.name)

            with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                    TimeElapsedColumn(),
                    console=console,
                    transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running DeepVariant for {sample_id}...", total=None)

                process = subprocess.Popen(shlex.split(f"su - {DOCKER_USER} -c '{deepvariant_command}'"),
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
                                           universal_newlines=True)

                for line in process.stdout:
                    progress.update(task, advance=1)
                    console.print(f"[dim]{line.strip()}[/dim]")
                    log_to_db(db, line.strip(), "DEBUG", "DeepVariant", sample_id, datadir.name)

                process.wait()

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(db, f"DeepVariant finished at {end_time}. Duration: {duration}", "INFO", "DeepVariant", sample_id,
                      datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"DeepVariant failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING",
                          "DeepVariant", sample_id, datadir.name)
            else:
                error_msg = f"DeepVariant failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "DeepVariant", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "DeepVariant", sample_id,
                          datadir.name)
                return "error"

        if output_vcf.exists() and output_gvcf.exists():
            console.print(Panel(f"[bold green]DeepVariant variants determined for {sample_id}[/bold green]"))
            log_to_api("DeepVariant variants determined", "INFO", "DeepVariant", sample_id, datadir.name)
            log_to_db(db, f"DeepVariant variants determined successfully for {sample_id}", "INFO", "DeepVariant",
                      sample_id, datadir.name)

            # Log file sizes and resource usage
            log_to_db(db, f"Output VCF size: {output_vcf.stat().st_size} bytes", "INFO", "DeepVariant", sample_id,
                      datadir.name)
            log_to_db(db, f"Output gVCF size: {output_gvcf.stat().st_size} bytes", "INFO", "DeepVariant", sample_id,
                      datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO",
                      "DeepVariant", sample_id, datadir.name)

            return "success"
        else:
            error_msg = f"DeepVariant completed but output VCF or gVCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "DeepVariant", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "DeepVariant", sample_id, datadir.name)
            return "error"

    return _deepvariant_caller()


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass