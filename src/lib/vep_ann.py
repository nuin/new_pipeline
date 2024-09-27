import subprocess
from pathlib import Path
from datetime import datetime
import psutil
import shlex
import shutil
import tempfile
import os

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.syntax import Syntax

from .log_api import log_to_api
from .db_logger import log_to_db, timer_with_db_log
from typing import Dict, Optional, Union
from tinydb import TinyDB

console = Console()

def get_vep_version(vep: str) -> str:
    """Get VEP version."""
    cmd = shlex.split(vep) + ["--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    for line in result.stdout.split('\n'):
        if "ensembl-vep" in line:
            return line.strip()
    return "Unknown"


def vep_annotate(
    sample_id: str,
    datadir: Path,
    vep: str,
    reference: Path,
    db: Union[TinyDB, Path],
    transcript_list: Optional[str] = None,
    max_retries: int = 3
) -> str:
    if not isinstance(db, TinyDB):
        raise TypeError(f"Expected TinyDB instance, got {type(db)}")

    @timer_with_db_log(db)
    def _annotate_merged():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        input_vcf = vcf_dir / f"{sample_id}_merged.vcf"
        output_vcf = vcf_dir / f"{sample_id}_merged.vep.vcf"
        vep_dir = Path("/apps/data/src/bin/vep")

        vep_version = get_vep_version(vep)
        log_to_db(db, f"Starting VEP annotation for sample {sample_id} with VEP version {vep_version}", "INFO", "vep_ann", sample_id, datadir.name)

        if not input_vcf.exists():
            error_msg = f"Input VCF file not found: {input_vcf}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
            return "error"

        if output_vcf.exists():
            console.print(Panel(f"[yellow]Annotated VCF file already exists for {sample_id}[/yellow]"))
            log_to_api("Annotated VCF file exists", "INFO", "vep_ann", sample_id, datadir.name)
            log_to_db(db, f"Annotated VCF file already exists for {sample_id}", "INFO", "vep_ann", sample_id, datadir.name)
            return "exists"

        console.print(Panel(f"[bold blue]Starting VEP annotation for {sample_id}[/bold blue]"))
        log_to_api("Start VEP annotation", "INFO", "vep_ann", sample_id, datadir.name)
        log_to_db(db, f"Starting VEP annotation for {sample_id}", "INFO", "vep_ann", sample_id, datadir.name)

        # Copy the input VCF to the VEP directory
        vep_input = vep_dir / input_vcf.name
        vep_output = vep_dir / output_vcf.name
        subprocess.run(['sudo', 'cp', str(input_vcf), str(vep_input)], check=True)

        # Change ownership and permissions of the input file
        subprocess.run(['sudo', 'chown', 'systemd-coredump:input', str(vep_input)], check=True)
        subprocess.run(['sudo', 'chmod', '644', str(vep_input)], check=True)

        vep_cmd = f"""
        sudo docker run --rm \
        -v {vep_dir}:/opt/vep/.vep:Z \
        ensemblorg/ensembl-vep \
        vep \
        -i /opt/vep/.vep/{input_vcf.name} \
        -o /opt/vep/.vep/{output_vcf.name} \
        --offline --cache --dir_cache /opt/vep/.vep \
        --assembly GRCh37 --species homo_sapiens \
        --fasta /opt/vep/.vep/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
        --format vcf --vcf --force_overwrite --hgvs \
        --dir /opt/vep/.vep --config /opt/vep/.vep/vep_config.ini
        """

        if transcript_list:
            vep_cmd = vep_cmd.replace("vep \\", f"vep --transcript_filter file=/opt/vep/.vep/{transcript_list} \\")

        console.print(Syntax(vep_cmd, "bash", theme="monokai", line_numbers=True))
        log_to_db(db, f"VEP command: {vep_cmd}", "INFO", "vep_ann", sample_id, datadir.name)

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(db, f"VEP started at {start_time} (Attempt {attempt + 1}/{max_retries})", "INFO", "vep_ann", sample_id, datadir.name)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console,
                transient=True
            ) as progress:
                task = progress.add_task(f"[cyan]Running VEP for {sample_id}...", total=None)

                process = subprocess.Popen(vep_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)

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
            log_to_db(db, f"VEP finished at {end_time}. Duration: {duration}", "INFO", "vep_ann", sample_id, datadir.name)

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(db, f"VEP failed. Retrying (Attempt {attempt + 2}/{max_retries})", "WARNING", "vep_ann", sample_id, datadir.name)
            else:
                error_msg = f"VEP failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
                log_to_db(db, f"{error_msg}. Last return code: {process.returncode}", "ERROR", "vep_ann", sample_id, datadir.name)
                return "error"

        if vep_output.exists():
            # Change ownership back to the original user
            subprocess.run(['sudo', 'chown', f'{os.getuid()}:{os.getgid()}', str(vep_output)], check=True)
            subprocess.run(['sudo', 'cp', str(vep_output), str(output_vcf)], check=True)
            subprocess.run(['sudo', 'rm', str(vep_input)], check=True)  # Remove the temporary input file
            subprocess.run(['sudo', 'rm', str(vep_output)], check=True)  # Remove the temporary output file
            console.print(Panel(f"[bold green]VEP annotation completed for {sample_id}[/bold green]"))
            log_to_api("VEP annotation completed", "INFO", "vep_ann", sample_id, datadir.name)
            log_to_db(db, f"VEP annotation completed successfully for {sample_id}", "INFO", "vep_ann", sample_id, datadir.name)

            log_to_db(db, f"Final VCF size: {output_vcf.stat().st_size} bytes", "INFO", "vep_ann", sample_id, datadir.name)
            log_to_db(db, f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB", "INFO", "vep_ann", sample_id, datadir.name)

            return "success"
        else:
            error_msg = f"VEP completed but output VCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "vep_ann", sample_id, datadir.name)
            return "error"

    return _annotate_merged()