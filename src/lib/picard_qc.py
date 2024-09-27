"""
.. module:: picard_qc
    :platform: Any
    :synopsis: Module that generates nucleotide-based coverage using Picard
.. moduleauthor:: Paulo Nuin, February 2015

"""

import shlex
import subprocess
from pathlib import Path
from typing import Union

from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from tinydb import TinyDB

from .db_logger import log_to_db
from .log_api import log_to_api

console = Console()


def get_coverage(
    sample_id: str,
    datadir: Path,
    reference: Path,
    picard: str,
    db: Union[TinyDB, str],
    bed_file: Union[str, Path],
    panel: str = "full",
) -> str:
    """
    .. module:: picard_qc
        :platform: Any
        :synopsis: Module that generates nucleotide-based coverage using Picard

    Function that generates nucleotide-based coverage using Picard CollectHsMetrics.

    :param sample_id: Sample identifier
    :type sample_id: str
    :param datadir: Path to the data directory
    :type datadir: Path
    :param reference: Path to the reference genome
    :type reference: Path
    :param picard: Path to the Picard jar file
    :type picard: str
    :param db: Database object or connection string
    :type db: Union[TinyDB, str]
    :param bed_file: Path to the BED file
    :type bed_file: Union[str, Path]
    :param panel: Type of panel (full or panel), defaults to "full"
    :type panel: str, optional
    :return: Status of the Picard CollectHsMetrics command
    :rtype: str
    """
    def _get_coverage():
        bam_file = datadir / "BAM" / sample_id / "BAM" / sample_id + ".bam"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)

        is_panel = panel != "full"
        output_suffix = ".nucl.panel.out" if is_panel else ".nucl.out"
        output_file = metrics_dir / f"{sample_id}{output_suffix}"
        metrics_output = metrics_dir / f"{sample_id}.{'panel.' if is_panel else ''}out"

        # Check if output file already exists and is not empty
        if output_file.exists() and output_file.stat().st_size > 0:
            message = f"Output file already exists: {output_file}"
            console.print(Panel(f"[bold yellow]{message}[/bold yellow]"))
            log_to_db(db, message, "INFO", "picard_coverage", sample_id, datadir.name)
            return "success"

        # Determine the correct intervals file
        if is_panel:
            intervals_file = Path(bed_file)
        else:
            # For full coverage, use the bait file
            intervals_file = Path(str(bed_file).replace(".bed", ".picard.bed"))

        picard_cmd = (
            f"java -Xmx8g -jar {picard} CollectHsMetrics "
            f"BI={intervals_file} "
            f"I={bam_file} "
            f"PER_BASE_COVERAGE={output_file} "
            f"MINIMUM_MAPPING_QUALITY=0 "
            f"MINIMUM_BASE_QUALITY=0 "
            f"TARGET_INTERVALS={intervals_file} "
            f"OUTPUT={metrics_output} "
            f"R={reference} "
            f"QUIET=true "
            f"TMP_DIR={metrics_dir} "
            f"USE_JDK_DEFLATER=true "
            f"USE_JDK_INFLATER=true "
            f"COMPRESSION_LEVEL=1 "
            f"MAX_RECORDS_IN_RAM=2000000 "
            f"VALIDATION_STRINGENCY=LENIENT"
        )

        console.print(Syntax(picard_cmd, "bash", theme="monokai", line_numbers=True))
        log_to_db(f"Picard command: {picard_cmd}", "INFO", "picard_coverage")

        try:
            process = subprocess.run(
                shlex.split(picard_cmd),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                universal_newlines=True,
            )

            for line in process.stdout.splitlines():
                console.print(f"[dim]{line.strip()}[/dim]")
                log_to_db(
                    db,
                    line.strip(),
                    "DEBUG",
                    "picard_coverage",
                    sample_id,
                    datadir.name,
                )

            if output_file.exists() and output_file.stat().st_size > 0:
                console.print(
                    Panel(
                        f"[bold green]Picard {panel} coverage file created for {sample_id}[/bold green]"
                    )
                )
                log_to_api(
                    f"Picard {panel} coverage file created",
                    "INFO",
                    "picard_coverage",
                    sample_id,
                    str(datadir),
                )
                log_to_db(
                    db,
                    f"Picard {panel} coverage file created successfully for {sample_id}",
                    "INFO",
                    "picard_coverage",
                    sample_id,
                    datadir.name,
                )
                return "success"
            else:
                raise FileNotFoundError(
                    f"Picard CollectHsMetrics completed but output file not found or empty for {sample_id}"
                )

        except subprocess.CalledProcessError as e:
            error_msg = f"Picard CollectHsMetrics failed for {sample_id}. Return code: {e.returncode}\nOutput: {e.output}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_coverage", sample_id, str(datadir))
            log_to_db(
                db, error_msg, "ERROR", "picard_coverage", sample_id, datadir.name
            )
            return "error"
        except Exception as e:
            error_msg = (
                f"Unexpected error in Picard CollectHsMetrics for {sample_id}: {str(e)}"
            )
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_coverage", sample_id, str(datadir))
            log_to_db(
                db, error_msg, "ERROR", "picard_coverage", sample_id, datadir.name
            )
            return "error"

    return _get_coverage()


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
