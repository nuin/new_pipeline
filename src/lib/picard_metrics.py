"""
.. module:: picard_metrics
    :platform: any
    :synopsis: Module that provides parsing capabilities for multiple types of Picard output
.. moduleauthor:: Paulo Nuin, January (modified May) 2016, Updated September 2024
"""

import subprocess
from pathlib import Path
from typing import Dict, Union

from rich.console import Console
from rich.markup import escape
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.syntax import Syntax
from tinydb import TinyDB

from .db_logger import log_to_db, timer_with_db_log
from .log_api import log_to_api

console = Console()


def run_picard_command(
    command: str,
    description: str,
    sample_id: str,
    datadir: Path,
    db: Union[TinyDB, str],
) -> int:
    """
    Runs a Picard command and logs output to a TinyDB database, and logs any errors to the API.

    Args:
        command: The Picard command to run.
        description: A description of the command being run.
        sample_id: The name of the sample being processed.
        datadir: The directory where data related to this sample is stored.
        db: A TinyDB database, or the path to a TinyDB database file.

    Returns:
        0 if the command runs successfully, 1 otherwise.
    """
    if isinstance(db, str):
        db = TinyDB(db)

    console.print(Panel(f"[bold blue]{description}[/bold blue]"))
    console.print(Syntax(command, "bash", theme="monokai", line_numbers=True))

    log_to_db(
        db,
        f"Running Picard command: {command}",
        "INFO",
        "picard_metrics",
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
            # Escape any potential markup in the output
            safe_line = escape(line.strip())
            console.print(f"[dim]{safe_line}[/dim]")
            log_to_db(
                db, line.strip(), "DEBUG", "picard_metrics", sample_id, datadir.name
            )

    process.wait()

    if process.returncode != 0:
        error_msg = f"Error in {description}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_api(error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
        log_to_db(db, error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
        return 1

    return 0


def get_yield(sample_id: str, datadir: Path, picard: str, db: Dict) -> str:
    """
    Runs Picard CollectQualityYieldMetrics on a sample and returns success or error.

    Args:
        sample_id: The name of the sample being processed.
        datadir: The directory where data related to this sample is stored.
        picard: The path to the Picard jar file.
        db: A dictionary representing the database.

    Returns:
        A string indicating whether the command was successful ("success"), an error occurred ("error"), or the file already exists ("exists").
    """

    @timer_with_db_log(db)
    def _get_yield():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        output_file = metrics_dir / f"{sample_id}.yield.out"

        if output_file.exists():
            console.print(
                Panel(
                    f"[yellow]Picard CollectQualityYieldMetrics file exists for {sample_id}[/yellow]"
                )
            )
            log_to_api(
                "Picard CollectQualityYieldMetrics file exists",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Picard CollectQualityYieldMetrics file exists for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "exists"

        picard_command = (
            f"{picard} CollectQualityYieldMetrics "
            f"I={bam_dir}/{sample_id}.bam "
            f"O={output_file}"
        )

        result = run_picard_command(
            picard_command,
            f"Generating Picard CollectQualityYieldMetrics for {sample_id}",
            sample_id,
            datadir,
            db,
        )

        if result == 0:
            console.print(
                Panel(
                    f"[bold green]CollectQualityYieldMetrics file created for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "CollectQualityYieldMetrics file created",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"CollectQualityYieldMetrics file created for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "success"
        else:
            return "error"

    return _get_yield()


def get_hs_metrics(
    sample_id: str,
    datadir: Path,
    reference: Path,
    bait_file: Union[str, Path],
    picard: str,
    db: Union[TinyDB, str],
    panel: str = "full",
) -> str:
    """
    Generate Picard CollectHsMetrics output for a given sample ID.

    Args:
    sample_id (str): The sample ID to generate the metrics for.
    datadir (Path): The base directory containing the sample's data.
    reference (Path): The path to the reference genome.
    bait_file (Union[str, Path]): The path to the bait file.
    picard (str): The path to the Picard jar file.
    db (Union[TinyDB, str]): The database to store the results in.
    panel (str): The panel to generate the metrics for. Defaults to "full".

    Returns:
    str: "exists" if the file already exists, "success" if the file was generated, or "error" if an error occurred.
    """

    @timer_with_db_log(db)
    def _get_hs_metrics():
        nonlocal bait_file, db  # Use the outer bait_file and db

        # Convert bait_file to Path if it's a string
        if isinstance(bait_file, str):
            bait_file = Path(bait_file)

        # If db is a string, create a TinyDB instance
        if isinstance(db, str):
            db = TinyDB(db)

        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"

        if panel == "full":
            output_file = metrics_dir / f"{sample_id}.hs_metrics.out"
        else:
            output_file = metrics_dir / f"{sample_id}.hs_metrics.panel.out"
            bait_file = bait_file.with_suffix(".picard.bed")

        if output_file.exists():
            console.print(
                Panel(
                    f"[yellow]Picard CollectHsMetrics {panel} file exists for {sample_id}[/yellow]"
                )
            )
            log_to_api(
                f"Picard CollectHsMetrics {panel} file exists",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Picard CollectHsMetrics {panel} file exists for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "exists"

        # Ensure bait_file exists
        if not bait_file.exists():
            error_msg = f"Bait file does not exist: {bait_file}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "picard_metrics", sample_id, datadir.name)
            return "error"

        picard_command = (
            f"{picard} CollectHsMetrics "
            f"-I {bam_dir}/{sample_id}.bam "
            f"-O {output_file} "
            f"-R {reference} "
            f"-BAIT_INTERVALS {bait_file} "
            f"-TARGET_INTERVALS {bait_file}"
        )

        result = run_picard_command(
            picard_command,
            f"Generating Picard CollectHsMetrics {panel} for {sample_id}",
            sample_id,
            datadir,
            db,
        )

        if result == 0:
            console.print(
                Panel(
                    f"[bold green]CollectHsMetrics {panel} file created for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                f"CollectHsMetrics {panel} file created",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"CollectHsMetrics {panel} file created for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "success"
        else:
            return "error"

    return _get_hs_metrics()


def get_align_summary(
    sample_id: str, datadir: Path, reference: Path, picard: str, db: Dict
) -> str:
    """
    Runs Picard's CollectAlignmentSummaryMetrics on a BAM file.

    If the output file already exists, logs a message and returns "exists".
    Otherwise, runs the command and logs a message with the result.
    Returns "success" or "error" depending on the result of the command.
    """

    @timer_with_db_log(db)
    def _get_align_summary():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        output_file = metrics_dir / f"{sample_id}.align_metrics.out"

        if output_file.exists():
            console.print(
                Panel(
                    f"[yellow]Picard AlignmentSummaryMetrics file exists for {sample_id}[/yellow]"
                )
            )
            log_to_api(
                "Picard AlignmentSummaryMetrics file exists",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Picard AlignmentSummaryMetrics file exists for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "exists"

        picard_command = (
            f"{picard} CollectAlignmentSummaryMetrics "
            f"I={bam_dir}/{sample_id}.bam "
            f"O={output_file} "
            f"R={reference}"
        )

        result = run_picard_command(
            picard_command,
            f"Generating Picard CollectAlignmentSummaryMetrics for {sample_id}",
            sample_id,
            datadir,
            db,
        )

        if result == 0:
            console.print(
                Panel(
                    f"[bold green]CollectAlignmentSummaryMetrics file created for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "CollectAlignmentSummaryMetrics file created",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"CollectAlignmentSummaryMetrics file created for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "success"
        else:
            return "error"

    return _get_align_summary()


def get_call_metrics(
    sample_id: str, datadir: Path, vcf_file: Path, picard: str, db: Dict
) -> str:
    """
    Run Picard CollectVariantCallingMetrics.

    This function runs the Picard CollectVariantCallingMetrics command to generate
    a variant calling metrics file for the given sample. If the file already exists,
    it will log a message and exit with "exists". If the command runs successfully,
    it will log a success message and return "success". If the command fails, it will
    log an error message and return "error".

    Parameters
    ----------
    sample_id : str
        The sample ID.
    datadir : Path
        The root directory for the pipeline output.
    vcf_file : Path
        The path to the DBSNP VCF file.
    picard : str
        The path to the Picard jar file.
    db : Dict
        The database connection object.

    Returns
    -------
    str
        The result of running the command. Either "exists", "success", or "error".
    """

    @timer_with_db_log(db)
    def _get_call_metrics():
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        output_file = metrics_dir / f"{sample_id}.call_metrics.out"
        input_vcf = vcf_dir / f"{sample_id}_merged.vcf"

        if (output_file.with_suffix(".variant_calling_detail_metrics")).exists():
            console.print(
                Panel(
                    f"[yellow]Picard CollectVariantCallingMetrics file exists for {sample_id}[/yellow]"
                )
            )
            log_to_api(
                "Picard CollectVariantCallingMetrics file exists",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"Picard CollectVariantCallingMetrics file exists for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "exists"

        picard_command = (
            f"{picard} CollectVariantCallingMetrics "
            f"I={input_vcf} "
            f"O={output_file} "
            f"DBSNP={vcf_file} "
            f"QUIET=true"
        )

        result = run_picard_command(
            picard_command,
            f"Generating Picard CollectVariantCallingMetrics for {sample_id}",
            sample_id,
            datadir,
            db,
        )

        if result == 0:
            console.print(
                Panel(
                    f"[bold green]CollectVariantCallingMetrics file created for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "CollectVariantCallingMetrics file created",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            log_to_db(
                db,
                f"CollectVariantCallingMetrics file created for {sample_id}",
                "INFO",
                "picard_metrics",
                sample_id,
                datadir.name,
            )
            return "success"
        else:
            return "error"

    return _get_call_metrics()


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass
