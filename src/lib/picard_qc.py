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
    bed_file: Path,
    panel: str = "full"
) -> str:
    def safe_log_to_db(message: str, level: str, program: str):
        if isinstance(db, TinyDB):
            log_to_db(db, message, level, program, sample_id, str(datadir))
        else:
            console.print(f"[yellow]Warning: DB logging unavailable. Message: {message}[/yellow]")

    def _get_coverage():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)

        is_panel = panel != "full"
        output_suffix = ".nucl.panel.out" if is_panel else ".nucl.out"
        output_file = metrics_dir / f"{sample_id}{output_suffix}"
        metrics_output = metrics_dir / f"{sample_id}.{'panel.' if is_panel else ''}out"
        bed_file_path = Path(bed_file) if isinstance(bed_file, str) else bed_file

        # Check if input files exist
        if not bam_dir.exists():
            error_msg = f"BAM directory not found: {bam_dir}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            safe_log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

        bam_file = bam_dir / f"{sample_id}.bam"
        if not bam_file.exists():
            error_msg = f"BAM file not found: {bam_file}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            safe_log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

        if not bed_file_path.exists():
            error_msg = f"BED file not found: {bed_file}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            safe_log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

        if not reference.exists():
            error_msg = f"Reference file not found: {reference}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            safe_log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"

        picard_cmd = (
            f"{picard} CollectHsMetrics "
            f"BI={bed_file_path} "
            f"I={bam_file} "
            f"PER_BASE_COVERAGE={output_file} "
            f"MINIMUM_MAPPING_QUALITY=0 "
            f"MINIMUM_BASE_QUALITY=0 "
            f"TARGET_INTERVALS={bed_file_path} "
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
        safe_log_to_db(f"Picard command: {picard_cmd}", "INFO", "picard_coverage")

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

            if output_file.exists():
                console.print(Panel(f"[bold green]Picard {panel} coverage file created for {sample_id}[/bold green]"))
                log_to_api(f"Picard {panel} coverage file created", "INFO", "picard_coverage", sample_id, str(datadir))
                safe_log_to_db(f"Picard {panel} coverage file created successfully for {sample_id}", "INFO", "picard_coverage")
                return "success"
            else:
                raise FileNotFoundError(f"Picard CollectHsMetrics completed but output file not found for {sample_id}")

        except subprocess.CalledProcessError as e:
            error_msg = f"Picard CollectHsMetrics failed for {sample_id}. Return code: {e.returncode}\nOutput: {e.output}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_coverage", sample_id, str(datadir))
            safe_log_to_db(error_msg, "ERROR", "picard_coverage")
            return "error"
        except Exception as e:
            error_msg = f"Unexpected error in Picard CollectHsMetrics for {sample_id}: {str(e)}"
            console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
            log_to_api(error_msg, "ERROR", "picard_coverage", sample_id, str(datadir))
            safe_log_to_db(error_msg, "ERROR", "picard_coverage")
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


# def chr_frame(segment: pd.DataFrame) -> tuple:
#     """
#     Function that checks chromosome regions for coverage under 25X
#
#     :param segment: current segment being analysed
#
#     :type segment: pandas DataFrame
#
#     :return: gene, segment and chromosome if there's a region under 25x, otherwise returns 'empty', 'empty', 'empty'
#
#     :rtype: tuple
#     """
#
#     if segment.loc[segment["coverage"].idxmin()]["coverage"] <= 100:
#         gene = segment.iloc[0]["Gene"]
#         chromosome = segment.iloc[0]["chrom"]
#         return gene, segment, chromosome
#
#     return "empty", "empty", "empty"


# def convert_g_to_c(chromosome, position, transcript):
#     """
#     Function that converts g. notation to c.
#
#     :param chromosome: chromosome number
#     :param position: position to be converted
#     :param transcript: list of transcripts
#
#     :type chromosome: string
#     :type position: integer
#     :type transcript: list
#
#     :return: HVGS notation of the position
#
#     """
#
#     fake_nucleotides = "A>A"
#     connection = Client(URL, cache=None)
#     mutalyzer = connection.service
#     result = mutalyzer.numberConversion(
#         build="hg19", variant=str(chromosome) + ":g." + str(position) + fake_nucleotides
#     )[0]
#     result_dict = {k: v for k, v in (x.split(":") for x in result)}
#     try:
#         return result_dict[transcript]
#     except Exception as e:
#         console.log(str(e))
#         return str(result[0])


# def create_table(under_30, chromosome, gene, trans, coverage_file):
#     """
#     Function that creates a text table with the locations under
#     20x of coverage
#
#     :param under30: list of regions under 100x (name should be changed)
#     :param chromosome: chromosome number
#     :param gene: gene symbol
#     :param trans: transcript ID
#     :param coverage_file: location of the coverage file
#
#     :type under30: list
#     :type chromosome: string
#     :type gene: string
#     :type transcript: string
#     :type coverage_file: string
#
#     :return: not return set
#
#     """
#
#     table_output = open(
#         os.path.dirname(coverage_file).replace("/Metrics", "")
#         + "/QC/"
#         + os.path.basename(coverage_file).replace(".nucl.out", "")
#         + "_under_25.txt",
#         "a",
#     )
#     for item in under_30:
#         if len(under_30[item]) != 0:
#             start_pos = convert_g_to_c(chromosome, under_30[item][0], trans)
#             end_pos = convert_g_to_c(chromosome, under_30[item][-1], trans)
#             table_output.write(
#                 chromosome
#                 + "\t"
#                 + gene
#                 + "\t"
#                 + str(int(under_30[item][-1]) - int(under_30[item][0]))
#                 + "\tg."
#             )
#             table_output.write(
#                 str(under_30[item][0]) + "\tg." + str(under_30[item][-1]) + "\t"
#             )
#             table_output.write(start_pos + "\t" + end_pos + "\n")
