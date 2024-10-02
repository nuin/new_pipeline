"""
.. module:: picard_qc
    :platform: Any
    :synopsis: Module that generates nucleotide-based coverage using Picard
.. moduleauthor:: Paulo Nuin, February 2015

"""

import subprocess
from pathlib import Path
from typing import Dict

from rich.console import Console
from rich.panel import Panel
from rich.progress import (BarColumn, Progress, SpinnerColumn, TextColumn,
                           TimeElapsedColumn)
from rich.syntax import Syntax

from .db_logger import log_to_db, timer_with_db_log
from .log_api import log_to_api

console = Console()


def run_picard_command(
    command: str, description: str, sample_id: str, datadir: Path, db: Dict
) -> int:
    console.print(Panel(f"[bold blue]{description}[/bold blue]"))
    console.print(Syntax(command, "bash", theme="monokai", line_numbers=True))

    log_to_db(
        db,
        f"Running Picard command: {command}",
        "INFO",
        "picard_qc",
        sample_id,
        datadir.name,
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task(description, total=None)

        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )

        for line in process.stdout:
            progress.update(task, advance=1)
            console.print(f"[dim]{line.strip()}[/dim]")
            log_to_db(db, line.strip(), "DEBUG", "picard_qc", sample_id, datadir.name)

    process.wait()

    if process.returncode != 0:
        error_msg = f"Error in {description}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_api(error_msg, "ERROR", "picard_qc", sample_id, datadir.name)
        log_to_db(db, error_msg, "ERROR", "picard_qc", sample_id, datadir.name)
        return 1

    return 0


def get_coverage(
    sample_id: str,
    datadir: Path,
    reference: Path,
    picard: str,
    db: Dict,
    bait_file: Path,
    panel: str = "full",
) -> str:
    @timer_with_db_log(db)
    def _get_coverage():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        bam_file = bam_dir / f"{sample_id}.bam"

        if panel == "full":
            output_file = metrics_dir / f"{sample_id}.coverage.out"
        else:
            output_file = metrics_dir / f"{sample_id}.coverage.panel.out"

        if output_file.exists():
            console.print(
                Panel(
                    f"[yellow]Picard coverage file already exists for {sample_id}[/yellow]"
                )
            )
            log_to_api(
                "Picard coverage file exists",
                "INFO",
                "picard_qc",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Picard coverage file exists for {sample_id}",
                "INFO",
                "picard_qc",
                sample_id,
                datadir.name,
            )
            return "exists"

        picard_command = (
            f"{picard} CollectWgsMetrics "
            f"I={bam_file} "
            f"O={output_file} "
            f"R={reference} "
            f"INTERVALS={bait_file} "
            f"COVERAGE_CAP=500 "
            f"INCLUDE_BQ_HISTOGRAM=true"
        )

        result = run_picard_command(
            picard_command,
            f"Generating Picard coverage metrics for {sample_id}",
            sample_id,
            datadir,
            db,
        )

        if result == 0:
            console.print(
                Panel(
                    f"[bold green]Picard coverage file created for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "Picard coverage file created",
                "INFO",
                "picard_qc",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Picard coverage file created for {sample_id}",
                "INFO",
                "picard_qc",
                sample_id,
                datadir.name,
            )
            return "success"
        else:
            error_msg = f"Failed to generate Picard coverage file for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_qc", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "picard_qc", sample_id, datadir.name)
            return "error"

    return _get_coverage()


# Other functions in picard_qc.py should be updated similarly


def get_transcripts(transcript_location: str) -> dict:
    """
    Function that gets the currently used transcripts for HGVS numberConversion

    :param transcript_location: Location of the transcript file

    :type transcript_location: string

    :return: Dictionary of transcripts

    :rtype: dict
    """

    transcript_file = open(transcript_location).read().splitlines()

    transcripts = {}
    for i in transcript_file:
        temp = i.split("\t")
        transcripts[temp[0]] = temp[1]

    return transcripts
