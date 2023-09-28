"""
.. module:: picard_qc
    :platform: Any
    :synopsis: Module that generates nucleotide-based coverage using Picard
.. moduleauthor:: Paulo Nuin, February 2015

"""

import os
import string
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd

from suds.client import Client

from pathlib import Path
from rich.console import Console

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
            return "exists"

        console.log(f"Starting Picard's CollectHsMetrics {sample_id}")
        picard_string = (
            f"{picard} CollectHsMetrics BI={bait_file} I={bam_dir}/{sample_id}.bam PER_BASE_COVERAGE={metrics_dir}/{sample_id}.nucl.out "
            f"MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 TARGET_INTERVALS={bait_file} OUTPUT={metrics_dir}/{sample_id}.out R={reference} QUIET=true"
        )
        proc = subprocess.Popen(
            picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        console.log(f"Picard command: {picard_string}")
        proc.wait()
        while True:
            output = proc.stderr.readline().strip()
            if output == b"":
                break
            else:
                console.log(output.decode("utf-8"))
        console.log(f"Picard coverage file created {sample_id}")
        return "success"

    if Path(f"{metrics_dir}/{sample_id}.nucl.panel.out").exists():
        console.log(f"Picard panel coverage file exists {sample_id}")
        return "exists"

    console.log(f"Starting Picard's CollectHsMetrics for panel {sample_id}")
    bait_file = bait_file.replace(".bed", ".picard.bed")
    picard_string = (
        f"{picard} CollectHsMetrics BI={bait_file} I={bam_dir}/{sample_id}.bam PER_BASE_COVERAGE={metrics_dir}/{sample_id}.nucl.panel.out "
        f"MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 TARGET_INTERVALS={bait_file} OUTPUT={metrics_dir}/{sample_id}.panel.out R={reference} QUIET=true"
    )
    console.log(f"Picard command: {picard_string}")
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
        logger.info("Picard coverage file exists")
        return "exists"

    logger.info("Starting Picard's CollectHsMetrics")
    picard_string = (
        '%s CollectHsMetrics BI="%s" I=%s.good.bam PER_BASE_COVERAGE=%s.nucl.out MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 TARGET_INTERVALS=%s OUTPUT=%s.out R=%s QUIET=true'
        % (picard, bait_file, argument2, argument, bait_file, argument, reference)
    )
    logger.info("Command " + picard_string)
    logger.info(picard_string)
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    while True:
        output = proc.stderr.readline().strip()
        if output == b"":
            break
        else:
            logger.info(output.decode("utf-8") + " " + sample_id)
    proc.wait()
    logger.info("INFO\t\tPicard coverage file created")
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
        logger.error(str(e) + " " + sample_id)
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


def generate_qc(coverage_file, sample_id, directory, transcript_location):
    """
    Function that manages all the work of the other functions
    calling each parser and generator and creating a pandas
    DataFrame for each consecutive location in the main coverage
    file

    :param coverage_file: location of the coverage file

    :type coverage_file: string

    """

    logger.info("Starting GenerateQC " + sample_id)

    transcripts = get_transcripts(transcript_location)
    coverage_info = pd.read_csv(
        coverage_file,
        sep="\t",
        header=(0),
        comment="*",
        usecols=range(4),
        dtype={"Total_Depth": int},
    )

    logger.info("Phase 1")
    try:
        empty, empty, empty, coverage_info["Gene"], empty = zip(
            *coverage_info["target"].apply(lambda x: x.split(":"))
        )
    except Exception as e:
        logger.warning(str(e))
        empty, coverage_info["Gene"], empty = zip(
            *coverage_info["target"].apply(lambda x: x.split(":"))
        )

    coverage_info["coverage"] = coverage_info["coverage"].astype(int)

    logger.info("Phase 2 - to plot")
    to_plot = defaultdict(list)
    for k, g in coverage_info.groupby(
        coverage_info["pos"] - np.arange(coverage_info.shape[0])
    ):
        gene, segment, chromosome = chr_frame(g)
        if gene != "empty":
            to_plot[gene].append(
                (
                    segment["coverage"].tolist(),
                    segment["pos"].tolist(),
                    segment["chrom"].tolist(),
                )
            )

    genes_plot(to_plot, coverage_file, transcripts, sample_id, directory)


if __name__ == "__main__":

    data_directory = "/Volumes/Jupiter/CancerPlusRuns/200604_NB551084_0096_AHCTWYAFX2_Cplus_2019_NGS_X2"
    sample_id = "17-310-020850_KS_OS"
    get_coverage(
        sample_id,
        data_directory,
        "/opt/reference/hg19.fasta",
        "/opt/BED/Inherited_Cancer_panel_BED_91122_Target_adjusted_FINAL_GIPoly.picard.bed",
        "picard",
        "/opt/bundle/transcripts.txt",
        "full",
    )
