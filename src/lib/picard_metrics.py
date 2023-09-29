"""
.. module:: picard_metrics
    :platform: any
    :synopsis: module that provides parsing capabilities for multiple types of picard output
.. moduleauthor:: Paulo Nuin, January (modified May) 2016

"""

import os
import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


def get_yield(sample_id, datadir, picard):
    """
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

    :todo: return error
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"
    metrics_dir = f"{datadir}/BAM/{sample_id}/Metrics/"

    if Path(f"{metrics_dir}{sample_id},yield.out").exists():
        console.log(f"Picard CollectQualityYieldMetrics file exists {sample_id}")
        return "exists"

    console.log(f"Starting Picard CollectQualityYieldMetrics creation {sample_id}")
    picard_string = (
        f"{picard} CollectQualityYieldMetrics I={bam_dir}{sample_id}.bam"
        f" O={metrics_dir}{sample_id}.yield.out"
    )
    console.log(f"{picard_string} {sample_id}")
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
    console.log("CollectQualityYieldMetrics file created " + sample_id)

    return "success"


def get_yield_parp(sample_id, directory, picard):
    """
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

    if os.path.isfile(argument + ".yield.out"):
        logger.info("Picard CollectQualityYieldMetrics file exists " + sample_id)
        return "exists"

    logger.info("Starting Picard CollectQualityYieldMetrics creation")
    picard_string = (
        "%s CollectQualityYieldMetrics I=%s.good.bam O=%s.yield.out QUIET=true"
        % (picard, argument2, argument)
    )
    logger.info(picard_string + " " + sample_id)
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    logger.info("CollectQualityYieldMetrics file created")
    return "success"


def get_hs_metrics(sample_id, directory, reference, bait_file, picard, panel="full"):
    """
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

    if panel == "full":
        if os.path.isfile(argument + ".hs_metrics.out") or os.path.isfile(
            argument + ".hs_metrics.out"
        ):
            logger.info("Picard CollectHsMetrics file exists " + sample_id)
            return "exists"

        logger.info("Generating Picard's CollectHsMetrics file " + sample_id)
        picard_string = (
            "%s CollectHsMetrics I=%s.recal_reads.bam O=%s.hs_metrics.out R=%s BAIT_INTERVALS=%s TARGET_INTERVALS=%s QUIET=true"
            % (picard, argument2, argument, reference, bait_file, bait_file)
        )
        proc = subprocess.Popen(
            picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        proc.wait()
        logger.info("CollectHsMetrics file created " + sample_id)
        return "success"

    bait_file = bait_file.replace(".bed", ".picard.bed")
    if os.path.isfile(argument + ".hs_metrics.panel.out") or os.path.isfile(
        argument2 + ".hs_metrics.panel.out"
    ):
        logger.info("Picard CollectHsMetrics panel file exists " + sample_id)
        return "exists"

    logger.info("Generating Picard CollectHsMetrics for panel " + sample_id)
    picard_string = (
        "%s CollectHsMetrics I=%s.recal_reads.bam O=%s.hs_metrics.panel.out R=%s BAIT_INTERVALS=%s TARGET_INTERVALS=%s QUIET=true"
        % (picard, argument2, argument, reference, bait_file, bait_file)
    )
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    logger.info(proc.stdout.read())
    logger.info("CollectHsMetrics file created " + sample_id)
    return "success"


def get_hs_metrics_parp(sample_id, directory, reference, bait_file, picard):
    """
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

    if os.path.isfile(argument + ".hs_metrics.out") or os.path.isfile(
        argument + ".hs_metrics.out"
    ):
        logger.info("Picard CollectHsMetrics file exists " + sample_id)
        return "exists"

    logger.info("Generating Picard's CollectHsMetrics file " + sample_id)
    picard_string = (
        "%s CollectHsMetrics I=%s.good.bam O=%s.hs_metrics.out R=%s BAIT_INTERVALS=%s TARGET_INTERVALS=%s QUIET=true"
        % (picard, argument2, argument, reference, bait_file, bait_file)
    )
    logger.info(picard_string + " " + sample_id)
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    logger.info("CollectHsMetrics file created " + sample_id)
    return "success"


def get_align_summary(sample_id, directory, reference, picard):
    """

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

    if os.path.isfile(argument + ".align_metrics.out"):
        logger.info("Picard CollectAlignmentSummaryMetrics file exists " + sample_id)
        return "exists"

    logger.info("Creating CollectAlignmentSummaryMetrics file " + sample_id)
    picard_string = (
        "%s CollectAlignmentSummaryMetrics I=%s.recal_reads.bam O=%s.align_metrics.out R=%s QUIET=true"
        % (picard, argument2, argument, reference)
    )
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    logger.info("CollectAlignmentSummaryMetrics file created " + sample_id)
    return "success"


def get_align_summary_parp(sample_id, directory, reference, picard):
    """

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

    if os.path.isfile(argument + ".align_metrics.out"):
        logger.info("Picard CollectAlignmentSummaryMetrics file exists " + sample_id)
        return "exists"

    logger.info("Creating CollectAlignmentSummaryMetrics file " + sample_id)
    picard_string = (
        "%s CollectAlignmentSummaryMetrics I=%s.good.bam O=%s.align_metrics.out R=%s QUIET=true"
        % (picard, argument2, argument, reference)
    )
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    logger.info("CollectAlignmentSummaryMetrics file created " + sample_id)
    return "success"


def get_call_metrics(sample_id, directory, vcf_file, picard):
    """
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

    :todo: return error
    :todo: fix argument
    """

    if os.path.isdir(directory + "/BAM/" + sample_id + "/Metrics/"):
        argument = directory + "/BAM/" + sample_id + "/Metrics/" + sample_id
    else:
        argument = directory + "/BAM/" + sample_id + "/" + sample_id

    if os.path.isdir(directory + "/BAM/" + sample_id + "/VCF/"):
        argument2 = directory + "/BAM/" + sample_id + "/VCF/" + sample_id
    else:
        argument2 = directory + "/BAM/" + sample_id + "/" + sample_id

    if os.path.isfile(argument + ".call_metrics.out.variant_calling_detail_metrics"):
        logger.info("Picard CollectVariantCallingMetrics file exists " + sample_id)
        return "exists"

    logger.info("Generating CollectVariantCallingMetrics file " + sample_id)
    picard_string = (
        "%s CollectVariantCallingMetrics I=%s_merged.vcf O=%s.call_metrics.out DBSNP=%s QUIET=true"
        % (picard, argument2, argument, vcf_file)
    )
    logger.info(picard_string + " " + sample_id)
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    logger.info("CollectVariantCallingMetrics file created " + sample_id)
    return "success"


if __name__ == "__main__":

    data_directory = "/Users/nuin/New_projects/Data/CAC0002-27826881"
    sample_id = "CAC0002"

    # get_yield(sample_id, data_directory, 'picard')
    # get_hs_metrics(sample_id, data_directory, '/opt/reference/hg19.fasta', '/opt/BED/Inherited_Cancer_panel_FINAL.picard.bed', '/opt/BED/Inherited_Cancer_panel_FINAL.list', 'picard')
    # get_align_summary(sample_id, data_directory, '/opt/reference/hg19.fasta', 'picard')
    # get_call_metrics(sample_id, data_directory, '/opt/bundle/dbsnp_138.hg19.vcf', 'picard')
