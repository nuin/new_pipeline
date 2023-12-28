"""
.. module:: picard_qc
    :platform: Any
    :synopsis: Module that generates nucleotide-based coverage using Picard
.. moduleauthor:: Paulo Nuin, February 2015

"""

import os
import string
import subprocess
from pathlib import Path

from rich.console import Console
from suds.client import Client

from .log_api import log_to_api

console = Console()
URL = "https://mutalyzer.nl/services/?wsdl"


class MyTemplate(string.Template):
    """
    TBA
    """

    delimiter = "%"
    idpattern = "[a-z]+_[a-z]+"


function = MyTemplate(
    """$(function () {
    $('#container%container_number').highcharts({
        chart:{
            type: 'spline'
        },
        title: {
            text: 'Coverage - %gene_name g.%location_here',

        },
        subtitle: {
            text: '',
        },
        xAxis: {
            categories: []
        },
        yAxis: {
            title: {
                text: 'coverage'
            },
            min: 0,
        },
        tooltip: {
            pointFormat: '{point.y}'
        },
        series: [{
            name: '%gene_name',
            data: %data_here,
            zones: [{
                    value: 26,
                    color: '#FF1811'
                    }],
        }]
    });
    });"""
)


def get_coverage(
    sample_id,
    datadir,
    reference,
    bait_file,
    picard,
    panel="full",
):
    """
    Function that calls Picard to generate nucleotide coverage

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

    :return: returns success or exists

    :todo: fix argument
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"
    metrics_dir = f"{datadir}/BAM/{sample_id}/Metrics/"

    if panel == "full":
        if Path(f"{metrics_dir}/{sample_id}.nucl.out").exists():
            console.log(f"Picard coverage file exists {sample_id}")
            log_to_api(
                "Picard coverage file exists",
                "INFO",
                "picard_coverage",
                sample_id,
                Path(datadir).name,
            )
            return "exists"

        console.log(f"Starting Picard's CollectHsMetrics {sample_id}")
        log_to_api(
            "Starting Picard's CollectHsMetrics",
            "INFO",
            "picard_coverage",
            sample_id,
            Path(datadir).name,
        )
        picard_string = (
            f"{picard} CollectHsMetrics BI={bait_file} I={bam_dir}/{sample_id}.bam PER_BASE_COVERAGE={metrics_dir}/{sample_id}.nucl.out "
            f"MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 TARGET_INTERVALS={bait_file} OUTPUT={metrics_dir}/{sample_id}.out R={reference} QUIET=true"
        )
        proc = subprocess.Popen(
            picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        console.log(f"Picard command: {picard_string}")
        log_to_api(
            picard_string, "INFO", "picard_coverage", sample_id, Path(datadir).name
        )
        proc.wait()
        while True:
            output = proc.stderr.readline().strip()
            if output == b"":
                break
            else:
                console.log(output.decode("utf-8"))
        console.log(f"Picard coverage file created {sample_id}")
        log_to_api(
            "Picard coverage file created",
            "INFO",
            "picard_coverage",
            sample_id,
            Path(datadir).name,
        )
        return "success"

    if Path(f"{metrics_dir}/{sample_id}.nucl.panel.out").exists():
        console.log(f"Picard panel coverage file exists {sample_id}")
        log_to_api(
            "Picard panel coverage file exists",
            "INFO",
            "picard_coverage",
            sample_id,
            Path(datadir).name,
        )
        return "exists"

    console.log(f"Starting Picard's CollectHsMetrics for panel {sample_id}")
    log_to_api(
        "Starting Picard's CollectHsMetrics for panel",
        "INFO",
        "picard_coverage",
        sample_id,
        Path(datadir).name,
    )
    bait_file = bait_file.replace(".bed", ".picard.bed")
    picard_string = (
        f"{picard} CollectHsMetrics BI={bait_file} I={bam_dir}/{sample_id}.bam PER_BASE_COVERAGE={metrics_dir}/{sample_id}.nucl.panel.out "
        f"MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 TARGET_INTERVALS={bait_file} OUTPUT={metrics_dir}/{sample_id}.panel.out R={reference} QUIET=true"
    )
    console.log(f"Picard command: {picard_string}")
    log_to_api(picard_string, "INFO", "picard_coverage", sample_id, Path(datadir).name)
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            console.log(output.decode("utf-8"))
    proc.wait()

    return "success"


def get_coverage_parp(sample_id, directory, reference, bait_file, picard):
    """
    Function that calls Picard to generate nucleotide coverage

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

    :return: returns success or exists
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


def get_transcripts(transcript_location):
    """
    Function that gets the currently used transcripts for HGVS numberConversion

    return: list of transcripts
    """

    transcript_file = open(transcript_location).read().splitlines()

    transcripts = {}
    for i in transcript_file:
        temp = i.split("\t")
        transcripts[temp[0]] = temp[1]

    return transcripts


def chr_frame(segment):
    """
    Function that checks chromosome regions for coverage under 25X

    :param segment: current segment being analysed
    :param transcripts: list of transcripts

    :type segment: pandas DataFrame
    :type transcripts: list

    :return: gene, segment and chromosome if there's a region under 25x

    """

    if segment.loc[segment["coverage"].idxmin()]["coverage"] <= 100:
        gene = segment.iloc[0]["Gene"]
        chromosome = segment.iloc[0]["chrom"]
        return gene, segment, chromosome

    return "empty", "empty", "empty"


def convert_g_to_c(chromosome, position, transcript):
    """
    Function that converts g. notation to c.

    :param chromosome: chromosome number
    :param position: position to be converted
    :param transcript: list of transcripts

    :type chromosome: string
    :type position: integer
    :type transcript: list

    :return: HVGS notation of the position

    """

    fake_nucleotides = "A>A"
    connection = Client(URL, cache=None)
    mutalyzer = connection.service
    result = mutalyzer.numberConversion(
        build="hg19", variant=str(chromosome) + ":g." + str(position) + fake_nucleotides
    )[0]
    result_dict = {k: v for k, v in (x.split(":") for x in result)}
    try:
        return result_dict[transcript]
    except Exception as e:
        console.log(str(e))
        return str(result[0])


def create_table(under_30, chromosome, gene, trans, coverage_file):
    """
    Function that creates a text table with the locations under
    20x of coverage

    :param under30: list of regions under 100x (name should be changed)
    :param chromosome: chromosome number
    :param gene: gene symbol
    :param trans: transcript ID
    :param coverage_file: location of the coverage file

    :type under30: list
    :type chromosome: string
    :type gene: string
    :type transcript: string
    :type coverage_file: string

    :return: not return set

    """

    table_output = open(
        os.path.dirname(coverage_file).replace("/Metrics", "")
        + "/QC/"
        + os.path.basename(coverage_file).replace(".nucl.out", "")
        + "_under_25.txt",
        "a",
    )
    for item in under_30:
        if len(under_30[item]) != 0:
            start_pos = convert_g_to_c(chromosome, under_30[item][0], trans)
            end_pos = convert_g_to_c(chromosome, under_30[item][-1], trans)
            table_output.write(
                chromosome
                + "\t"
                + gene
                + "\t"
                + str(int(under_30[item][-1]) - int(under_30[item][0]))
                + "\tg."
            )
            table_output.write(
                str(under_30[item][0]) + "\tg." + str(under_30[item][-1]) + "\t"
            )
            table_output.write(start_pos + "\t" + end_pos + "\n")
