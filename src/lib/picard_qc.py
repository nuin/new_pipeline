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


# class MyTemplate(string.Template):
#     """
#     TBA
#     """
#
#     delimiter = "%"
#     idpattern = "[a-z]+_[a-z]+"
#
#
# function = MyTemplate(
#     """$(function () {
#     $('#container%container_number').highcharts({
#         chart:{
#             type: 'spline'
#         },
#         title: {
#             text: 'Coverage - %gene_name g.%location_here',
#
#         },
#         subtitle: {
#             text: '',
#         },
#         xAxis: {
#             categories: []
#         },
#         yAxis: {
#             title: {
#                 text: 'coverage'
#             },
#             min: 0,
#         },
#         tooltip: {
#             pointFormat: '{point.y}'
#         },
#         series: [{
#             name: '%gene_name',
#             data: %data_here,
#             zones: [{
#                     value: 26,
#                     color: '#FF1811'
#                     }],
#         }]
#     });
#     });"""
# )




def get_picard_version(picard: str) -> str:
    """Get Picard version."""
    result = subprocess.run(["java", "-jar", picard, "--version"], capture_output=True, text=True)
    return result.stdout.strip()


def get_coverage(sample_id: str, datadir: Path, reference: Path, bait_file: Path, picard: str, db: Dict,
                 panel: str = "full", threads: int = 4, max_retries: int = 3) -> str:
    @timer_with_db_log(db)
    def _get_coverage():
        bam_dir = datadir / "BAM" / sample_id / "BAM"
        metrics_dir = datadir / "BAM" / sample_id / "Metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)

        is_panel = panel != "full"
        output_suffix = ".nucl.panel.out" if is_panel else ".nucl.out"
        output_file = metrics_dir / f"{sample_id}{output_suffix}"
        metrics_output = metrics_dir / f"{sample_id}.{'panel.' if is_panel else ''}out"

        picard_version = get_picard_version(picard)
        log_to_db(db, f"Starting Picard CollectHsMetrics for sample {sample_id} with Picard version {picard_version}",
                  "INFO", "picard_coverage", sample_id, datadir.name)

        if output_file.exists():
            console.print(Panel(f"[yellow]Picard coverage file already exists for {sample_id}[/yellow]"))
            log_to_api(f"Picard {panel} coverage file exists", "INFO", "picard_coverage", sample_id, datadir.name)
            log_to_db(db, f"Picard {panel} coverage file already exists for {sample_id}", "INFO", "picard_coverage",
                      sample_id, datadir.name)
            return "exists"

        console.print(
            Panel(f"[bold blue]Starting Picard's CollectHsMetrics for {panel} coverage of {sample_id}[/bold blue]"))
        log_to_api(f"Starting Picard's CollectHsMetrics for {panel} coverage", "INFO", "picard_coverage", sample_id,
                   datadir.name)
        log_to_db(db, f"Starting Picard's CollectHsMetrics for {panel} coverage of {sample_id}", "INFO",
                  "picard_coverage", sample_id, datadir.name)

        # Ensure bait_file is a Path object
        bait_file = Path(bait_file)

        # Modify bait_file for panel if necessary
        if is_panel:
            bait_file = bait_file.with_suffix('.picard.bed')

        picard_cmd = (
            f"java -jar {picard} CollectHsMetrics "
            f"BI={bait_file} "
            f"I={bam_dir}/{sample_id}.bam "
            f"PER_BASE_COVERAGE={output_file} "
            f"MINIMUM_MAPPING_QUALITY=0 "
            f"MINIMUM_BASE_QUALITY=0 "
            f"TARGET_INTERVALS={bait_file} "
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
        log_to_db(db, f"Picard command: {picard_cmd}", "INFO", "picard_coverage", sample_id, datadir.name)

        # ... (rest of the function remains the same)

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
