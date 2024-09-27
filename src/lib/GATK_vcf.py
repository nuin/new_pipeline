"""
.. module:: GATK_vcf
    :platform: any
    :synopsis: Module that calls GATK to compare VCF files and performs post-analysis of these variants
.. moduleauthor:: Paulo Nuin, April 2016, Updated September 2024
"""

import shlex
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict

import psutil
from rich.console import Console
from rich.panel import Panel
from rich.progress import (BarColumn, Progress, SpinnerColumn, TextColumn,
                           TimeElapsedColumn)
from rich.syntax import Syntax

from .db_logger import log_to_db, timer_with_db_log
from .log_api import log_to_api

console = Console()


def get_gatk_version(gatk: str) -> str:
    """
    Retrieves the version of the Genome Analysis Toolkit (GATK) executable.

    Parameters:
        gatk (str): The path to the GATK executable.

    Returns:
        str: The version of the GATK executable.
    """
    result = subprocess.run(
        ["java", "-jar", gatk, "--version"], capture_output=True, text=True
    )
    return result.stdout.strip()


def filter_merged_vcf(
    input_vcf: Path, output_vcf: Path, db: Dict, sample_id: str, datadir: Path
) -> bool:
    """
    Filters a merged VCF file based on the number of callers and the presence of the Intersection set.

    Parameters:
        input_vcf (Path): The path to the input VCF file.
        output_vcf (Path): The path to the output VCF file.
        db (Dict): A dictionary containing database connection information.
        sample_id (str): The ID of the sample being processed.
        datadir (Path): The path to the directory containing the data.

    Returns:
        bool: True if the filtering was successful, False otherwise.
    """
    try:
        with input_vcf.open("r") as infile, output_vcf.open("w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    outfile.write(line)
                    continue

                fields = line.strip().split("\t")
                info = fields[7]

                if "set=" not in info:
                    continue

                set_field = [f for f in info.split(";") if f.startswith("set=")][0]
                callers = set_field.split("=")[1].split("-")

                # Keep variant if called by at least two callers
                # or if it's in the Intersection set
                if len(callers) >= 2 or "Intersection" in callers:
                    outfile.write(line)
                else:
                    log_to_db(
                        db,
                        f"Filtered out variant: {fields[0]}:{fields[1]} {fields[3]}>{fields[4]} (set={'-'.join(callers)})",
                        "DEBUG",
                        "GATK_vcf",
                        sample_id,
                        datadir.name,
                    )

        log_to_db(
            db,
            f"Successfully filtered merged VCF for {sample_id}",
            "INFO",
            "GATK_vcf",
            sample_id,
            datadir.name,
        )
        return True
    except Exception as e:
        error_msg = f"Error filtering merged VCF for {sample_id}: {str(e)}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_api(error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
        log_to_db(db, error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
        return False


def vcf_comparison(
    datadir: Path,
    sample_id: str,
    reference: Path,
    gatk: str,
    db: Dict,
    max_retries: int = 3,
) -> str:
    """
    Performs a comparison of VCF files using GATK, merging and filtering the results.

    Parameters:
        datadir (Path): The directory containing the data.
        sample_id (str): The ID of the sample being processed.
        reference (Path): The reference genome file.
        gatk (str): The path to the GATK executable.
        db (Dict): The database connection.
        max_retries (int): The maximum number of retries for the GATK command. Defaults to 3.

    Returns:
        str: The status of the operation, either "exists", "success", or "error".
    """
    @timer_with_db_log(db)
    def _vcf_comparison():
        vcf_dir = datadir / "BAM" / sample_id / "VCF"
        final_vcf = vcf_dir / f"{sample_id}_merged.vcf"

        gatk_version = get_gatk_version(gatk)
        log_to_db(
            db,
            f"Starting VCF comparison for sample {sample_id} with GATK version {gatk_version}",
            "INFO",
            "GATK_vcf",
            sample_id,
            datadir.name,
        )

        if final_vcf.exists():
            console.print(
                Panel(
                    f"[yellow]Merged VCF file already exists for {sample_id}[/yellow]"
                )
            )
            log_to_api(
                "Merged VCF file exists", "INFO", "GATK_vcf", sample_id, datadir.name
            )
            log_to_db(
                db,
                f"Merged VCF file already exists for {sample_id}",
                "INFO",
                "GATK_vcf",
                sample_id,
                datadir.name,
            )
            return "exists"

        console.print(
            Panel(
                f"[bold blue]Starting merge of GATK, GATK3, Freebayes, and Octopus VCFs for {sample_id}[/bold blue]"
            )
        )
        log_to_api(
            "Starting merge of GATK, GATK3, Freebayes, and Octopus VCFs",
            "INFO",
            "GATK_vcf",
            sample_id,
            datadir.name,
        )
        log_to_db(
            db,
            f"Starting merge of GATK, GATK3, Freebayes, and Octopus VCFs for {sample_id}",
            "INFO",
            "GATK_vcf",
            sample_id,
            datadir.name,
        )

        gatk_command = (
            f"{gatk} -T CombineVariants "
            f"-R {reference} "
            f"--variant:freebayes {vcf_dir}/{sample_id}_freebayes.vcf "
            f"--variant:gatk {vcf_dir}/{sample_id}_GATK.vcf "
            f"--variant:gatk3 {vcf_dir}/{sample_id}_GATK3.vcf "
            f"--variant:octopus {vcf_dir}/{sample_id}_octopus.vcf "
            f"-o {final_vcf} "
            f"--genotypemergeoption UNSORTED "
            f"--mergeInfoWithMaxAC "
            f"--minimumN 2"
        )

        console.print(Syntax(gatk_command, "bash", theme="monokai", line_numbers=True))
        log_to_db(
            db,
            f"GATK command: {gatk_command}",
            "INFO",
            "GATK_vcf",
            sample_id,
            datadir.name,
        )

        for attempt in range(max_retries):
            start_time = datetime.now()
            log_to_db(
                db,
                f"GATK CombineVariants started at {start_time} (Attempt {attempt + 1}/{max_retries})",
                "INFO",
                "GATK_vcf",
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
                task = progress.add_task(
                    f"[cyan]Running GATK CombineVariants for {sample_id}...", total=None
                )

                process = subprocess.Popen(
                    shlex.split(gatk_command),
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
                        db, line.strip(), "DEBUG", "GATK_vcf", sample_id, datadir.name
                    )

                process.wait()

            end_time = datetime.now()
            duration = end_time - start_time
            log_to_db(
                db,
                f"GATK CombineVariants finished at {end_time}. Duration: {duration}",
                "INFO",
                "GATK_vcf",
                sample_id,
                datadir.name,
            )

            if process.returncode == 0:
                break
            elif attempt < max_retries - 1:
                log_to_db(
                    db,
                    f"GATK CombineVariants failed. Retrying (Attempt {attempt + 2}/{max_retries})",
                    "WARNING",
                    "GATK_vcf",
                    sample_id,
                    datadir.name,
                )
            else:
                error_msg = f"GATK CombineVariants failed for {sample_id} after {max_retries} attempts"
                console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
                log_to_api(error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
                log_to_db(
                    db,
                    f"{error_msg}. Last return code: {process.returncode}",
                    "ERROR",
                    "GATK_vcf",
                    sample_id,
                    datadir.name,
                )
                return "error"

        if final_vcf.exists():
            console.print(
                Panel(
                    f"[bold green]VCFs merged successfully for {sample_id}[/bold green]"
                )
            )
            log_to_api(
                "VCFs merged successfully", "INFO", "GATK_vcf", sample_id, datadir.name
            )
            log_to_db(
                db,
                f"VCFs merged successfully for {sample_id}",
                "INFO",
                "GATK_vcf",
                sample_id,
                datadir.name,
            )

            # Filter the merged VCF
            console.print(
                Panel(f"[bold blue]Filtering merged VCF for {sample_id}[/bold blue]")
            )
            log_to_api(
                "Filtering merged VCF", "INFO", "GATK_vcf", sample_id, datadir.name
            )
            log_to_db(
                db,
                f"Filtering merged VCF for {sample_id}",
                "INFO",
                "GATK_vcf",
                sample_id,
                datadir.name,
            )

            with tempfile.NamedTemporaryFile(
                mode="w+t", delete=False, suffix=".vcf", dir=vcf_dir
            ) as temp_file:
                temp_vcf = Path(temp_file.name)
                if filter_merged_vcf(final_vcf, temp_vcf, db, sample_id, datadir):
                    # Replace the original merged VCF with the filtered one
                    shutil.move(str(temp_vcf), str(final_vcf))
                    console.print(
                        Panel(
                            f"[bold green]VCF filtered and replaced successfully for {sample_id}[/bold green]"
                        )
                    )
                    log_to_api(
                        "VCF filtered and replaced successfully",
                        "INFO",
                        "GATK_vcf",
                        sample_id,
                        datadir.name,
                    )
                    log_to_db(
                        db,
                        f"VCF filtered and replaced successfully for {sample_id}",
                        "INFO",
                        "GATK_vcf",
                        sample_id,
                        datadir.name,
                    )

                    # Log file size and resource usage
                    log_to_db(
                        db,
                        f"Final VCF size: {final_vcf.stat().st_size} bytes",
                        "INFO",
                        "GATK_vcf",
                        sample_id,
                        datadir.name,
                    )
                    log_to_db(
                        db,
                        f"Peak memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB",
                        "INFO",
                        "GATK_vcf",
                        sample_id,
                        datadir.name,
                    )

                    return "success"
                else:
                    return "error"
        else:
            error_msg = f"GATK CombineVariants completed but output VCF not found for {sample_id}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
            log_to_db(db, error_msg, "ERROR", "GATK_vcf", sample_id, datadir.name)
            return "error"

    return _vcf_comparison()


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass
