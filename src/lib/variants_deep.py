import logging
import shlex
import subprocess
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Union

import psutil
from rich.console import Console
from rich.panel import Panel
from rich.progress import (BarColumn, Progress, SpinnerColumn, TextColumn,
                           TimeElapsedColumn)
from rich.syntax import Syntax

from .db_logger import log_to_db, timer_with_db_log
from .log_api import log_to_api

console = Console()

# Set up a fallback logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DEEPVARIANT_IMAGE = "google/deepvariant:1.6.1"
TIMEOUT = 3600  # 1 hour timeout


def safe_log_to_db(
    db: Optional[Dict],
    message: str,
    level: str,
    program: str,
    sample_id: str,
    run_id: str,
):
    """Safely log to database, falling back to console logging if db is None."""
    if db is not None:
        log_to_db(db, message, level, program, sample_id, run_id)
    else:
        logger.info(f"{level} - {program} - {sample_id} - {run_id} - {message}")


def run_with_sudo(
    command: str, timeout: int = TIMEOUT, sample_id: str = "", db: Optional[Dict] = None
) -> subprocess.CompletedProcess:
    """Run a command with sudo and timeout, capturing ongoing output."""
    full_command = f"sudo {command}"

    def log_output(pipe, log_func):
        for line in iter(pipe.readline, ""):
            log_func(line.strip())

    try:
        process = subprocess.Popen(
            shlex.split(full_command),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )

        stdout_thread = threading.Thread(
            target=log_output,
            args=(
                process.stdout,
                lambda msg: safe_log_to_db(
                    db, msg, "INFO", "DeepVariant", sample_id, ""
                ),
            ),
        )
        stderr_thread = threading.Thread(
            target=log_output,
            args=(
                process.stderr,
                lambda msg: safe_log_to_db(
                    db, msg, "ERROR", "DeepVariant", sample_id, ""
                ),
            ),
        )

        stdout_thread.start()
        stderr_thread.start()

        process.wait(timeout=timeout)

        stdout_thread.join()
        stderr_thread.join()

        return subprocess.CompletedProcess(
            args=full_command, returncode=process.returncode, stdout="", stderr=""
        )
    except subprocess.TimeoutExpired:
        console.print(
            Panel(
                f"[bold red]Command timed out after {timeout} seconds: {command}[/bold red]"
            )
        )
        return subprocess.CompletedProcess(
            args=full_command, returncode=-1, stdout="", stderr="Timeout"
        )


def deepvariant_caller(
    datadir: Path,
    sample_id: str,
    reference: Union[str, Path],
    bed_file: Union[str, Path],
    db: Optional[Dict],
    threads: int = 16,
    max_retries: int = 3,
) -> str:
    @timer_with_db_log(db)
    def _deepvariant_caller():
        # Convert string paths to Path objects
        reference_path = Path(reference) if isinstance(reference, str) else reference
        bed_file_path = Path(bed_file) if isinstance(bed_file, str) else bed_file

        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_vcf = vcf_dir / f"{sample_id}_DeepVariant.vcf.gz"
        output_gvcf = vcf_dir / f"{sample_id}_DeepVariant.g.vcf.gz"
        input_bam = bam_dir / f"{sample_id}.bam"
        intermediate_results_dir = vcf_dir / "deepvariant_intermediate"

        safe_log_to_db(
            db,
            f"Starting DeepVariant for sample {sample_id} with DeepVariant",
            "INFO",
            "DeepVariant",
            sample_id,
            datadir.name,
        )

        if output_vcf.exists() and output_gvcf.exists():
            console.print(
                Panel(
                    f"[yellow]DeepVariant VCF and gVCF files already exist for {sample_id}[/yellow]"
                )
            )
            log_to_api(
                "DeepVariant VCF and gVCF files exist",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )
            safe_log_to_db(
                db,
                f"DeepVariant VCF and gVCF files already exist for {sample_id}",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )
            return "exists"

        console.print(
            Panel(
                f"[bold blue]Starting variant calling with DeepVariant for {sample_id}[/bold blue]"
            )
        )
        log_to_api(
            "Start variant calling with DeepVariant",
            "INFO",
            "DeepVariant",
            sample_id,
            datadir.name,
        )
        safe_log_to_db(
            db,
            f"Starting variant calling with DeepVariant for {sample_id}",
            "INFO",
            "DeepVariant",
            sample_id,
            datadir.name,
        )

        # Ensure output directories exist
        vcf_dir.mkdir(parents=True, exist_ok=True)
        intermediate_results_dir.mkdir(parents=True, exist_ok=True)

        # Prepare Docker volume mounts
        volume_mounts = [
            f"-v {datadir}:/input",
            f"-v {vcf_dir}:/output",
            f"-v {reference_path.parent}:/reference",
            f"-v {bed_file_path.parent}:/regions",
        ]

        deepvariant_command = (
            f"docker run --rm {' '.join(volume_mounts)} {DEEPVARIANT_IMAGE} "
            f"/opt/deepvariant/bin/run_deepvariant "
            f"--model_type=WGS "
            f"--ref=/reference/{reference_path.name} "
            f"--reads=/input/{input_bam.relative_to(datadir)} "
            f"--regions=/regions/{bed_file_path.name} "
            f"--output_vcf=/output/{output_vcf.name} "
            f"--output_gvcf=/output/{output_gvcf.name} "
            f"--intermediate_results_dir=/output/{intermediate_results_dir.name} "
            f"--num_shards={threads} "
            f"--logging_level=debug "
            f"--verbose"
        )

        console.print(
            Syntax(deepvariant_command, "bash", theme="monokai", line_numbers=True)
        )
        safe_log_to_db(
            db,
            f"DeepVariant command: {deepvariant_command}",
            "INFO",
            "DeepVariant",
            sample_id,
            datadir.name,
        )

        for attempt in range(max_retries):
            start_time = datetime.now()
            safe_log_to_db(
                db,
                f"DeepVariant started at {start_time} (Attempt {attempt + 1}/{max_retries})",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )

            result = run_with_sudo(
                deepvariant_command, timeout=TIMEOUT, sample_id=sample_id, db=db
            )

            end_time = datetime.now()
            duration = end_time - start_time
            safe_log_to_db(
                db,
                f"DeepVariant finished at {end_time}. Duration: {duration}",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )

            if result.returncode == 0:
                break
            elif result.returncode == -1:
                error_msg = (
                    f"DeepVariant timed out for {sample_id} after {TIMEOUT} seconds"
                )
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "DeepVariant", sample_id, datadir.name)
                safe_log_to_db(
                    db, error_msg, "ERROR", "DeepVariant", sample_id, datadir.name
                )
                return "error"
            elif attempt < max_retries - 1:
                safe_log_to_db(
                    db,
                    f"DeepVariant failed. Retrying (Attempt {attempt + 2}/{max_retries})",
                    "WARNING",
                    "DeepVariant",
                    sample_id,
                    datadir.name,
                )
            else:
                error_msg = (
                    f"DeepVariant failed for {sample_id} after {max_retries} attempts"
                )
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "DeepVariant", sample_id, datadir.name)
                safe_log_to_db(
                    db,
                    f"{error_msg}. Last return code: {result.returncode}",
                    "ERROR",
                    "DeepVariant",
                    sample_id,
                    datadir.name,
                )
                console.print(
                    Panel(f"[bold red]Error output: {result.stderr}[/bold red]")
                )
                return "error"

        if output_vcf.exists() and output_gvcf.exists():
            console.print(
                Panel(
                    f"[bold green]DeepVariant variants determined for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "DeepVariant variants determined",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )
            safe_log_to_db(
                db,
                f"DeepVariant variants determined successfully for {sample_id}",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )

            # Log file sizes and resource usage
            safe_log_to_db(
                db,
                f"Output VCF size: {output_vcf.stat().st_size} bytes",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )
            safe_log_to_db(
                db,
                f"Output gVCF size: {output_gvcf.stat().st_size} bytes",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )
            safe_log_to_db(
                db,
                f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB",
                "INFO",
                "DeepVariant",
                sample_id,
                datadir.name,
            )

            return "success"
        else:
            error_msg = f"DeepVariant completed but output VCF or gVCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "DeepVariant", sample_id, datadir.name)
            safe_log_to_db(
                db, error_msg, "ERROR", "DeepVariant", sample_id, datadir.name
            )
            return "error"

    return _deepvariant_caller()


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass
