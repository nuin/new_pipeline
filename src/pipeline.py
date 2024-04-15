"""
.. module:: align
    :platform: any
    :synopsis: This module is the main piece in the NGS pipeline, performing
                FASTQ alignments and calling all other steps in the process
.. moduleauthor:: Paulo Nuin, December 2015

"""

import glob
import os
from collections import defaultdict
from pathlib import Path

import click
import yaml
from dotenv import dotenv_values
from rich.console import Console

from lib.bwa_align import run_bwa
from lib.cnv import cnv_calculation, compile_samples
from lib.count2 import extract_counts
from lib.dup_indels import remove_duplicates
from lib.enrichment import get_enrichment
from lib.extract_identity import create_identity_table, mpileup
from lib.GATK_vcf import vcf_comparison
from lib.log_api import log_to_api
from lib.picard_actions import picard_sort
from lib.picard_metrics import get_align_summary, get_hs_metrics, get_yield
from lib.picard_qc import get_coverage
from lib.process_identity import barcoding, compile_barcodes
from lib.recalibration import base_recal1, recalibrate
from lib.snpEff_ann import annotate_merged
from lib.uniformity import get_coverage_values
from lib.utils import compile_identity, move_bam
from lib.variants_freebayes import edit_freebayes_vcf, freebayes_caller
from lib.variants_GATK import haplotype_caller
from lib.variants_GATK3 import haplotype_caller as haplotype_caller3
from lib.variants_octopus import octopus_caller

# main configuration file
# couch_credentials = open('lib/config/couchdb').read().splitlines()
# run configuration file - should be moved to some argument parsing lib

# checks current version
try:
    VERSIONFILE = "VERSION"
    VERSION = open(VERSIONFILE).read().strip()
except:
    VERSION = "N/A" 

console = Console()


def split_string(ctx, param, value):
    """ """

    if value is None:
        return []
    else:
        return [item.strip() for item in value.split(",")]


def find_fastq(datadir, panel_samples, panel):
    """

    :param datadir:
    :return:
    """

    fastqs = glob.glob(f"{datadir}/*.fastq.gz")
    fastqs = sorted(fastqs)

    to_return = []
    for fastq in fastqs:
        console.log(f"Found {fastq}")
        log_to_api(f"Found {fastq}", "INFO", "pipeline", "NA", Path(datadir).name)
        to_return.append(Path(fastq).name)

    return fastqs


def get_ids(fastqs):
    """

    :param fastqs:
    :return:
    """

    sample_ids = []
    for fastq in fastqs:
        sample_ids.append(Path(fastq).name.split("_")[0])

    sample_ids = sorted(set(sample_ids))

    return sample_ids


def create_directories(datadir, sample_ids, panel):
    """

    :param datadir:
    :param fastqs:
    :return:
    """

    try:
        os.makedirs(f"{datadir}/BAM/")
    except FileExistsError:
        console.log(f"{datadir}/BAM/ already exists")
    console.log(f"Creating directories in {panel} directory")
    for sample in sample_ids:
        console.log(f"Processing {sample}")
        log_to_api(
            f"Processing {sample}", "INFO", "pipeline", sample, Path(datadir).name
        )
        if not Path(f"{datadir}/BAM/{sample}").exists():
            console.log(f"Creating {datadir}/BAM/{sample}")
            log_to_api(
                f"Creating {datadir}/BAM/{sample}",
                "INFO",
                "pipeline",
                sample,
                Path(datadir).name,
            )
            os.makedirs(f"{datadir}/BAM/{sample}")
        if not Path(f"{datadir}/BAM/{sample}/BAM/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/BAM/")
            log_to_api(
                f"Creating {datadir}/BAM/{sample}/BAM/",
                "INFO",
                "pipeline",
                sample,
                Path(datadir).name,
            )
            os.makedirs(f"{datadir}/BAM/{sample}/BAM/")
        if not Path(f"{datadir}/BAM/{sample}/VCF/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/VCF/")
            log_to_api(
                f"Creating {datadir}/BAM/{sample}/VCF/",
                "INFO",
                "pipeline",
                sample,
                Path(datadir).name,
            )
            os.makedirs(f"{datadir}/BAM/{sample}/VCF/")
        if not Path(f"{datadir}/BAM/{sample}/QC/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/QC/")
            log_to_api(
                f"Creating {datadir}/BAM/{sample}/QC/",
                "INFO",
                "pipeline",
                sample,
                Path(datadir).name,
            )
            os.makedirs(f"{datadir}/BAM/{sample}/QC/")
        if not Path(f"{datadir}/BAM/{sample}/Metrics/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/Metrics/")
            log_to_api(
                f"Creating {datadir}/BAM/{sample}/Metrics/",
                "INFO",
                "pipeline",
                sample,
                Path(datadir).name,
            )
            os.makedirs(f"{datadir}/BAM/{sample}/Metrics/")

    console.log("Directories created/existed")


def align_files(config, datadir: Path, samples: list, fastqs: list) -> bool:
    """

    :param datadir:
    :param samples:
    :return:
    """

    env = dotenv_values(f"{Path.cwd()}/.env")
    configuration = yaml.safe_load(open(config))
    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    console.log("Aligning files")
    log_to_api("Aligning files", "INFO", "pipeline", "NA", Path(datadir).name)

    for sample in samples:
        console.log(f"Processing {sample}")
        log_to_api(
            f"Processing {sample}", "INFO", "pipeline", sample, Path(datadir).name
        )
        fastq_files = glob.glob(f"{datadir}/{sample}*.fastq.gz")
        console.log(f"Aligning sample {sample}")
        log_to_api(
            f"Aligning sample {sample}", "INFO", "pipeline", sample, Path(datadir).name
        )
        run_bwa(sample, fastq_files, datadir, reference, bwa, samtools)


def process_dir(config, datadir, samples, panel):
    """
    Function that runs all the steps of the pipeline, usually calling
    functions from other modules

    :param datadir: Datadir with the location of the FASTQ files

    :type datadir: string

    :return: to_return a dictionary with the fail/success of each step
    :rtype: dictionary

    :todo: decrease the number of arguments
    :todo: replace Picard coverage with GATK
    """

    env = dotenv_values(f"{Path.cwd()}/.env")
    configuration = yaml.safe_load(open(config))

    console.log("Looking for FASTQ files")
    log_to_api("Looking for FASTQ files", "INFO", "pipeline", "NA", Path(datadir).name)
    panel_samples = list(configuration["BED"].keys())
    fastqs = find_fastq(datadir, panel_samples, panel)
    sample_ids = get_ids(fastqs)

    console.log(f"Found {len(sample_ids)} samples")
    log_to_api(
        f"Found {len(sample_ids)} samples", "INFO", "pipeline", "NA", Path(datadir).name
    )
    console.log(f"Found {' '.join(sample_ids)}")
    log_to_api(
        f"Found {' '.join(sample_ids)}", "INFO", "pipeline", "NA", Path(datadir).name
    )

    if len(samples) >= 1:
        console.log(f"Samples to be analysed: {' '.join(samples)}")
        log_to_api(
            f"Samples to be analysed: {' '.join(samples)}",
            "INFO",
            "pipeline",
            "NA",
            Path(datadir).name,
        )
        sample_ids = samples

    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    console.log(f"Using BWA: {bwa}")
    log_to_api(f"Using BWA: {bwa}", "INFO", "pipeline", "NA", Path(datadir).name)
    console.log(f"Using SAMTOOLS: {samtools}")
    log_to_api(
        f"Using SAMTOOLS: {samtools}", "INFO", "pipeline", "NA", Path(datadir).name
    )
    console.log(f"Using reference: {reference}")
    log_to_api(
        f"Using reference: {reference}", "INFO", "pipeline", "NA", Path(datadir).name
    )

    create_directories(datadir, sample_ids, panel)
    align_files(config, datadir, sample_ids, fastqs)
    analyse_pairs(config, datadir, sample_ids, panel)

    return True


def analyse_pairs(config, datadir, samples, panel):
    """
    Main function of the file. Gets all pairs (samples) to
    be analysed and guide the whole process, step by step

    :param config yaml file for the run
    :param datadir location of the run files
    :param software_conf yaml file that list the software to be used
    :param pairs list of samples in the run
    :param yaml_file run configuration file (required for CNV, needs checking)


    :type config: string
    :type datadir: string
    :type software_conf: string
    :type pairs: dict
    :type yaml_file: string

    :todo: better return
    :todo full refactoring
    :todo check yaml_file paramneter
    """

    to_return = defaultdict(dict)
    env = dotenv_values(f"{Path.cwd()}/.env")
    configuration = yaml.safe_load(open(config))

    # load software configuration
    samtools = env["SAMTOOLS"]
    picard = env["PICARD"]
    gatk = env["GATK"]
    gatk3 = env["GATK3"]
    freebayes = env["FREEBAYES"]
    snpEff = env["SNPEFF"]
    octopus = env["OCTOPUS"]
    reference = configuration["Reference"]
    vcf_file = configuration["VCF"]
    bed_file = configuration["BED"]
    bait_file = configuration["BAIT"]

    for pos, sample in enumerate(samples):
        console.log(f"Processing {sample} :: {pos + 1} of {len(samples)}")
        log_to_api(
            f"Processing {sample} :: {pos + 1} of {len(samples)}",
            "INFO",
            "pipeline",
            sample,
            Path(datadir).name,
        )
        rm_duplicates = remove_duplicates(sample, datadir, picard)
        to_return[sample]["dedup"] = rm_duplicates
        move_bam(datadir, sample, "dedup")

        recalibration_step1 = base_recal1(
            datadir, sample, bed_file[sample], vcf_file, reference, gatk
        )
        to_return[sample]["recalibration1"] = recalibration_step1

        recalibration_final = recalibrate(datadir, sample, reference, gatk)
        to_return[sample]["recalibrate"] = recalibration_final
        move_bam(datadir, sample, "recal_reads")

        haplotype_caller(datadir, sample, reference, bed_file[sample], gatk)

        haplotype_caller3(datadir, sample, reference, bed_file[sample], gatk3)

        freebayes_caller(datadir, sample, reference, bed_file[sample], freebayes)
        to_return[sample]["picard_sort"] = picard_sort(
            datadir, sample, reference, picard
        )

        to_return[sample]["freebayes_edit"] = edit_freebayes_vcf(sample, datadir)

        to_return[sample]["variants_octopus"] = octopus_caller(
            datadir, sample, reference, bed_file[sample], octopus
        )

        to_return[sample]["vcf_merge"] = vcf_comparison(
            datadir, sample, reference, gatk3
        )

        to_return[sample]["snpEff"] = annotate_merged(sample, datadir, snpEff)

        to_return[sample]["picard_coverage"] = get_coverage(
            sample, datadir, reference, bait_file, picard
        )

        to_return[sample]["picard_coverage_panel"] = get_coverage(
            sample, datadir, reference, bed_file[sample], picard, "panel"
        )

        to_return[sample]["picard_yield"] = get_yield(sample, datadir, picard)

        to_return[sample]["picard_hs_metrics"] = get_hs_metrics(
            sample, datadir, reference, bait_file, picard
        )
        to_return[sample]["picard_hs_metrics_panel"] = get_hs_metrics(
            sample, datadir, reference, bed_file[sample], picard, "panel"
        )

        to_return[sample]["picard_align_metrics"] = get_align_summary(
            sample, datadir, reference, picard
        )

        to_return[sample]["mpileup_ident"] = mpileup(
            sample, datadir, "/apps/data/src/bundle/identity.txt", samtools
        )
        to_return[sample]["identity_table"] = create_identity_table(sample, datadir)
        to_return[sample]["full_identity"] = barcoding(sample, datadir)

        if panel == "Cplus":
            to_return[sample]["cnv"] = extract_counts(
                datadir,
                "/apps/data/src/BED/new/C+_ALL_IDPE_01JUN2021_Window.bed",
                sample,
            )
        else:
            to_return[sample]["cnv"] = extract_counts(
                datadir,
                "/apps/data/src/BED/new/CardiacALL_29MAR2021_Window.bed",
                sample,
            )

        console.log(f"Processing {sample} completed")
        log_to_api(
            f"Processing {sample} completed",
            "INFO",
            "pipeline",
            sample,
            Path(datadir).name,
        )

    console.log("Compiling identity file")
    log_to_api("Compiling identity file", "INFO", "pipeline", "NA", Path(datadir).name)
    if not Path(f"{datadir}/identity.txt").exists():
        console.log("Identity file does not exist, creating it")
        log_to_api(
            "Identity file does not exist, creating it",
            "INFO",
            "pipeline",
            "NA",
            Path(datadir).name,
        )
        compile_identity(datadir)
    else:
        log_to_api("Identity file exists", "INFO", "pipeline", "NA", Path(datadir).name)
        console.log("Identity file exists")

    console.log("Compiling barcodes")
    log_to_api("Compiling barcodes", "INFO", "pipeline", "NA", Path(datadir).name)
    compile_barcodes(datadir)

    console.log("Calculating and saving CNV read normalization")
    log_to_api(
        "Calculating and saving CNV read normalization",
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )
    all_cnvs = compile_samples(datadir)
    cnv_calculation(datadir, all_cnvs, config)

    console.log("Calculating uniformity")
    log_to_api("Calculating uniformity", "INFO", "pipeline", "NA", Path(datadir).name)
    get_coverage_values(datadir, panel)

    console.log("Calculating enrichment")
    log_to_api("Calculating enrichment", "INFO", "pipeline", "NA", Path(datadir).name)
    for pos, sample in enumerate(samples):
        get_enrichment(sample, datadir, panel)

    return to_return


def generate_analysis(config, datadir, samples, panel):
    """
    Main function of the script, find input files, start alignemt and processing

    :param datadir: location of the FASTQ files

    :type datadir: string

    :return: patients dictionary
    :rtype: dictionary
    """

    console.log("Starting analysis and qc report generation")
    log_to_api(
        "Starting analysis and qc report generation",
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )
    console.log(f"Configuration file: {config}")
    log_to_api(
        f"Configuration file: {config}", "INFO", "pipeline", "NA", Path(datadir).name
    )
    console.log(f"Datadir: {datadir}")
    log_to_api(f"Datadir: {datadir}", "INFO", "pipeline", "NA", Path(datadir).name)
    console.log(f"Panel: {panel}")
    log_to_api(f"Panel: {panel}", "INFO", "pipeline", "NA", Path(datadir).name)
    if samples == []:
        console.log("Samples: ALL")
        log_to_api("Samples: ALL", "INFO", "pipeline", "NA", Path(datadir).name)
    else:
        console.log(f"Samples: {' '.join(samples)}")
        log_to_api(
            f"Samples: {' '.join(samples)}",
            "INFO",
            "pipeline",
            "NA",
            Path(datadir).name,
        )

    s = process_dir(config, datadir, samples, panel)

    return s


@click.command()
@click.option(
    "-c", "--configuration-file", "configuration_file", help="YAML file", required=True
)
@click.option(
    "-s",
    "--sample",
    "samples",
    help="sample to be analysed (it can be multiple, each with a -s)",
    callback=split_string,
    required=False,
)
@click.option("-d", "--datadir", "datadir", help="run directory", required=True)
@click.option("-p", "--panel", "panel", help="panel to be used", required=True)
def run_analysis(configuration_file, datadir, panel, samples):
    """

    :param configuration_file:
    :param samples:
    :param datadir:
    :return:
    """
    if not samples:
        console.log("No samples provided, will analyse all samples in the run")
        samples = []
    console.log("Pipeline current version is " + VERSION)
    log_to_api(
        "Pipeline current version is " + VERSION,
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )
    console.log("All requirements found, starting analysis")
    log_to_api(
        "All requirements found, starting analysis",
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )

    sample_dict = generate_analysis(configuration_file, datadir, samples, panel)

    finished = open(f"/nfs/mgn_dna/NGS/temp/automation/finished_{panel}", "w")
    finished.write(f"{datadir}\n")
    finished.close()

    return sample_dict


if __name__ == "__main__":
    run_analysis()

# 23GN-135G00063
# 23GN-135G00072
# 23GN-125G00029
