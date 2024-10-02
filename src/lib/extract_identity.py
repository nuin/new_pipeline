"""
.. module:: extract_identity
    :platform: Any
    :synopsis: This module extracts information on the nucleotide in the 16 locations that determine identity
.. moduleauthor:: Paulo Nuin, August 2016, Updated September 2024
"""

import subprocess
from collections import Counter
from pathlib import Path
from typing import Dict, Tuple

from rich.console import Console
from rich.panel import Panel
from rich.progress import (BarColumn, Progress, SpinnerColumn, TextColumn,
                           TimeElapsedColumn)
from rich.syntax import Syntax

from .db_logger import log_to_db, timer_with_db_log
from .log_api import log_to_api

console = Console()


def run_command(
    command: str, description: str, sample_id: str, datadir: Path, db: Dict
) -> int:
    """
    Runs a shell command and logs output to the database.

    Args:
        command (str): The shell command to run
        description (str): A description of the command to display in the progress bar
        sample_id (str): The sample ID
        datadir (Path): The path to the directory containing the sample's data
        db (Dict): The database connection

    Returns:
        int: 0 if the command was successful, 1 if there was an error
    """
    console.print(Panel(f"[bold blue]{description}[/bold blue]"))
    console.print(Syntax(command, "bash", theme="monokai", line_numbers=True))

    log_to_db(
        db,
        f"Running command: {command}",
        "INFO",
        "extract_identity",
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
            log_to_db(
                db, line.strip(), "DEBUG", "extract_identity", sample_id, datadir.name
            )

    process.wait()

    if process.returncode != 0:
        error_msg = f"Error in {description}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_api(error_msg, "ERROR", "extract_identity", sample_id, datadir.name)
        log_to_db(db, error_msg, "ERROR", "extract_identity", sample_id, datadir.name)
        return 1

    return 0


def mpileup(
    sample_id: str, datadir: Path, identity: Path, samtools: str, db: Dict
) -> str:
    """
    Generate a mpileup file for the given sample_id and datadir using the given samtools executable and identity file.

    Args:
        sample_id (str): The sample ID
        datadir (Path): The directory containing the sample files
        identity (Path): The path to the identity file
        samtools (str): The path to the samtools executable
        db (Dict): The database connection

    Returns:
        str: One of 'exists', 'success', or 'error'
    """

    @timer_with_db_log(db)
    def _mpileup():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        output_file = bam_dir / "identity.mpileup"
        input_bam = bam_dir / f"{sample_id}.bam"

        if output_file.exists():
            console.print(
                Panel(f"[yellow]Identity mpileup file exists for {sample_id}[/yellow]")
            )
            log_to_api(
                "Identity mpileup file exists",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Identity mpileup file exists for {sample_id}",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            return "exists"

        console.print(
            Panel(
                f"[bold blue]Starting mpileup process for identity file {sample_id}[/bold blue]"
            )
        )
        log_to_api(
            "Starting mpileup process for identity file",
            "INFO",
            "extract_identity",
            sample_id,
            datadir.name,
        )
        log_to_db(
            db,
            f"Starting mpileup process for identity file {sample_id}",
            "INFO",
            "extract_identity",
            sample_id,
            datadir.name,
        )

        mpileup_command = (
            f"{samtools} mpileup -l {identity} {input_bam} > {output_file}"
        )

        result = run_command(
            mpileup_command,
            f"Generating mpileup for identity of {sample_id}",
            sample_id,
            datadir,
            db,
        )

        if result == 0:
            console.print(
                Panel(
                    f"[bold green]mpileup generation completed for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "mpileup generation completed",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"mpileup generation completed for {sample_id}",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            return "success"
        else:
            return "error"

    return _mpileup()


def get_nucleotides(pileup: str) -> Tuple[int, int, int, int]:
    """
    Counts the number of each nucleotide in a pileup string.

    Args:
        pileup (str): The pileup string to count

    Returns:
        Tuple[int, int, int, int]: A tuple of (A, C, G, T) counts
    """
    counter = Counter(char.upper() for char in pileup if char.isalpha())
    return counter["A"], counter["C"], counter["G"], counter["T"]


def create_identity_table(sample_id: str, datadir: Path, db: Dict) -> str:
    """
    Creates a tab-delimited identity file from the mpileup file generated by `create_mpileup`.

    The identity file has the following columns:
        - chromosome
        - position
        - A count
        - C count
        - G count
        - T count

    The function returns "exists" if the identity file already exists, "success" if the file is created successfully, and "error" if there is an error.

    :param sample_id: The sample ID of the BAM file.
    :param datadir: The directory containing the BAM file.
    :param db: The database connection.
    :return: A string indicating whether the file was created successfully or not.
    """

    @timer_with_db_log(db)
    def _create_identity_table():
        mpileup_file = datadir / "BAM" / sample_id / "BAM" / "identity.mpileup"
        identity_file = datadir / "BAM" / sample_id / "identity.txt"

        if identity_file.exists():
            console.print(
                Panel(f"[yellow]Identity file exists for {sample_id}[/yellow]")
            )
            log_to_api(
                "Identity file exists",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Identity file exists for {sample_id}",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            return "exists"

        console.print(
            Panel(f"[bold blue]Creating identity file for {sample_id}[/bold blue]")
        )
        log_to_api(
            "Creating identity file",
            "INFO",
            "extract_identity",
            sample_id,
            datadir.name,
        )
        log_to_db(
            db,
            f"Creating identity file for {sample_id}",
            "INFO",
            "extract_identity",
            sample_id,
            datadir.name,
        )

        try:
            with mpileup_file.open("r") as f, identity_file.open("w") as identity:
                for line in f:
                    chrom, pos, _, _, pileup, _ = line.split("\t")
                    nucleotides = get_nucleotides(pileup)
                    identity.write(
                        f"{chrom}\t{pos}\t{nucleotides[0]}\t{nucleotides[1]}\t{nucleotides[2]}\t{nucleotides[3]}\n"
                    )

            console.print(
                Panel(
                    f"[bold green]Identity file creation complete for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "Identity file creation complete",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Identity file creation complete for {sample_id}",
                "INFO",
                "extract_identity",
                sample_id,
                datadir.name,
            )
            return "success"
        except Exception as e:
            error_msg = f"Error creating identity file for {sample_id}: {str(e)}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "extract_identity", sample_id, datadir.name)
            log_to_db(
                db, error_msg, "ERROR", "extract_identity", sample_id, datadir.name
            )
            return "error"

    return _create_identity_table()


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass
