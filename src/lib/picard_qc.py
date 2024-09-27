"""
.. module:: picard_qc
    :platform: Any
    :synopsis: Module that generates nucleotide-based coverage using Picard
.. moduleauthor:: Paulo Nuin, February 2015

"""

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


URL = "https://mutalyzer.nl/services/?wsdl"

import shlex
import subprocess
from pathlib import Path
from typing import Union
from tinydb import TinyDB

def get_coverage(
    sample_id: str,
    datadir: Path,
    reference: Path,
    picard: str,
    db: Union[TinyDB, str],
    bed_file: Union[str, Path],
    panel: str = "full"
) -> str:
    def _get_coverage():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
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
            log_to_db(message, "INFO", "picard_coverage")
            return "success"

        # Determine the correct intervals file
        if is_panel:
            intervals_file = Path(bed_file)
        else:
            # For full coverage, use the bait file
            intervals_file = Path(str(bed_file).replace('.bed', '.picard.bed'))

        # Check if input files exist
        if not bam_dir.exists():
            error_msg = f"BAM directory not found: {bam_dir}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

        bam_file = bam_dir / f"{sample_id}.bam"
        if not bam_file.exists():
            error_msg = f"BAM file not found: {bam_file}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            safe_log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

        if not intervals_file.exists():
            error_msg = f"Intervals file not found: {intervals_file}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

        if not reference.exists():
            error_msg = f"Reference file not found: {reference}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

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
                universal_newlines=True
            )

            for line in process.stdout.splitlines():
                console.print(f"[dim]{line.strip()}[/dim]")
                safe_log_to_db(line.strip(), "DEBUG", "picard_coverage")

            if output_file.exists() and output_file.stat().st_size > 0:
                console.print(Panel(f"[bold green]Picard {panel} coverage file created for {sample_id}[/bold green]"))
                log_to_api(f"Picard {panel} coverage file created", "INFO", "picard_coverage", sample_id, str(datadir))
                safe_log_to_db(f"Picard {panel} coverage file created successfully for {sample_id}", "INFO", "picard_coverage")
                return "success"
            else:
                raise FileNotFoundError(f"Picard CollectHsMetrics completed but output file not found or empty for {sample_id}")

        except subprocess.CalledProcessError as e:
            error_msg = f"Picard CollectHsMetrics failed for {sample_id}. Return code: {e.returncode}\nOutput: {e.output}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_coverage", sample_id, str(datadir))
            log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"
        except Exception as e:
            error_msg = f"Unexpected error in Picard CollectHsMetrics for {sample_id}: {str(e)}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_coverage", sample_id, str(datadir))
            log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

    return _get_coverage()



def get_coverage_parp(
    sample_id: str, directory: str, reference: str, bait_file: str, picard: str
) -> str:
    """
    Function that calls Picard to generate nucleotide coverage for PARP.

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param directory: Location of the BAM files
    :param reference: Reference genome
    :param bait_file: Picard specific BED file
    :param picard: Picard jar file location

    :type sample_id: string
    :type directory: string
    :type reference: string
    :type bait_file: string
    :type picard: string

    :return: returns 'success' if the Picard coverage file is successfully created, 'exists' if the file already exists.

    :todo: return error
    """

    if os.path.isdir(directory + "/BAM/" + sample_id + "/Metrics/"):
        argument = directory + "/BAM/" + sample_id + "/Metrics/" + sample_id
    else:
        argument = directory + "/BAM/" + sample_id + "/" + sample_id

    if os.path.isdir(directory + "/BAM/" + sample_id + "/BAM/"):
        argument2 = directory + "/BAM/" + sample_id + "/BAM/" + sample_id
    else:
        argument2 = directory + "/BAM/" + sample_id + "/" + sample_id

    if os.path.isfile(argument + ".nucl.out"):
        return "exists"

    picard_string = (
        '%s CollectHsMetrics BI="%s" I=%s.good.bam PER_BASE_COVERAGE=%s.nucl.out MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 TARGET_INTERVALS=%s OUTPUT=%s.out R=%s QUIET=true'
        % (picard, bait_file, argument2, argument, bait_file, argument, reference)
    )
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            pass
    proc.wait()
    return "success"


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
