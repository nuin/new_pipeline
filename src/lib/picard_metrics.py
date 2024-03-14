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

from .log_api import log_to_api

console = Console()


def get_yield(sample_id: str, datadir: str, picard: str) -> str:
    """
    Function that generates Picard CollectQualityYieldMetrics.

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param datadir: Location of the BAM files
    :param picard: Picard jar file location

    :type sample_id: string
    :type datadir: string
    :type picard: string

    :return: returns 'success' if the Picard CollectQualityYieldMetrics file is successfully created, 'exists' if the file already exists.

    :todo: return error
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"
    metrics_dir = f"{datadir}/BAM/{sample_id}/Metrics/"

    if Path(f"{metrics_dir}{sample_id}.yield.out").exists():
        console.log(f"Picard CollectQualityYieldMetrics file exists {sample_id}")
        log_to_api
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


def get_hs_metrics(
    sample_id: str,
    datadir: str,
    reference: str,
    bait_file: str,
    picard: str,
    panel: str = "full",
) -> str:
    """
    Function that generates Picard CollectHsMetrics.

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param datadir: Location of the BAM files
    :param reference: Reference genome
    :param bait_file: Picard specific BED file
    :param picard: Picard jar file location
    :param panel: Panel type, default is "full"

    :type sample_id: string
    :type datadir: string
    :type reference: string
    :type bait_file: string
    :type picard: string
    :type panel: string

    :return: returns 'success' if the Picard CollectHsMetrics file is successfully created, 'exists' if the file already exists.

    :todo: return error
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"
    metrics_dir = f"{datadir}/BAM/{sample_id}/Metrics/"

    if panel == "full":
        if Path(f"{metrics_dir}{sample_id}.hs_metrics.out").exists():
            console.log(f"Picard CollectHsMetrics file exists {sample_id}")
            return "exists"

        console.log(f"Generating Picard CollectHsMetrics {sample_id}")

        picard_string = (
            f"{picard} CollectHsMetrics I={bam_dir}{sample_id}.bam"
            f" O={metrics_dir}{sample_id}.hs_metrics.out R={reference}"
            f" BAIT_INTERVALS={bait_file} TARGET_INTERVALS={bait_file}"
        )

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
        console.log("CollectHsMetrics file created " + sample_id)

        return "success"

    bait_file = bait_file.replace(".bed", ".picard.bed")
    if Path(f"{metrics_dir}{sample_id}.hs_metrics.panel.out").exists():
        console.log(f"Picard CollectHsMetrics panel file exists {sample_id}")
        return "exists"

    console.log(f"Generating Picard CollectHsMetrics panel {sample_id}")
    picard_string = (
        f"{picard} CollectHsMetrics I={bam_dir}{sample_id}.bam"
        f" O={metrics_dir}{sample_id}.hs_metrics.panel.out R={reference}"
        f" BAIT_INTERVALS={bait_file} TARGET_INTERVALS={bait_file}"
    )

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
    console.log("CollectHsMetrics panel file created " + sample_id)
    return "success"


def get_align_summary(sample_id: str, datadir: str, reference: str, picard: str) -> str:
    """
    Function that generates Picard CollectAlignmentSummaryMetrics.

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param datadir: Location of the BAM files
    :param reference: Reference genome
    :param picard: Picard jar file location

    :type sample_id: string
    :type datadir: string
    :type reference: string
    :type picard: string

    :return: returns 'success' if the Picard CollectAlignmentSummaryMetrics file is successfully created, 'exists' if the file already exists.

    :todo: return error
    """

    bam_dir = f"{datadir}/BAM/{sample_id}/BAM/"
    metrics_dir = f"{datadir}/BAM/{sample_id}/Metrics/"

    if Path(f"{metrics_dir}{sample_id}.align_metrics.out").exists():
        console.log(f"Picard  AlignSummary file exists {sample_id}")
        log_to_api(
            "Picard  AlignSummary file exists",
            "INFO",
            "picard_align_summary",
            sample_id,
            Path(datadir).name,
        )
        return "exists"

    console.log(f"Generating Picard CollectAlignmentSummaryMetrics {sample_id}")
    log_to_api(
        "Generating Picard CollectAlignmentSummaryMetrics",
        "INFO",
        "picard_align_summary",
        sample_id,
        Path(datadir).name,
    )
    picard_string = (
        f"{picard} CollectAlignmentSummaryMetrics I={bam_dir}{sample_id}.bam"
        f" O={metrics_dir}{sample_id}.align_metrics.out R={reference}"
    )
    console.log(f"{picard_string} {sample_id}")
    log_to_api(
        picard_string, "INFO", "picard_align_summary", sample_id, Path(datadir).name
    )
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    console.log(f"CollectAlignmentSummaryMetrics file created {sample_id}")
    log_to_api(
        "CollectAlignmentSummaryMetrics file created",
        "INFO",
        "picard_align_summary",
        sample_id,
        Path(datadir).name,
    )
    return "success"


def get_call_metrics(sample_id: str, directory: str, vcf_file: str, picard: str) -> str:
    """
    Function that generates Picard CollectVariantCallingMetrics.

    :param sample_id: ID of the patient/sample being analysed using Picard
    :param directory: Location of the BAM files
    :param vcf_file: VCF file used for variant calling
    :param picard: Picard jar file location

    :type sample_id: string
    :type directory: string
    :type vcf_file: string
    :type picard: string

    :return: returns 'success' if the Picard CollectVariantCallingMetrics file is successfully created, 'exists' if the file already exists.

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
        return "exists"

    picard_string = (
        "%s CollectVariantCallingMetrics I=%s_merged.vcf O=%s.call_metrics.out DBSNP=%s QUIET=true"
        % (picard, argument2, argument, vcf_file)
    )
    proc = subprocess.Popen(
        picard_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.wait()
    return "success"
