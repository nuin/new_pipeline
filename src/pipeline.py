"""
.. module:: align
    :platform: any
    :synopsis: This module is the main piece in the NGS pipeline, performing
                FASTQ alignments and calling all other steps in the process
.. moduleauthor:: Paulo Nuin, December 2015

"""

import glob
import logging
import os
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import click
import yaml
from rich.console import Console

console = Console()

from dotenv import dotenv_values

from lib.bwa_align import run_bwa
from lib.dup_indels import remove_duplicates
from lib.recalibration import base_recal1, recalibrate
from lib.variants_GATK import haplotype_caller
from lib.variants_GATK3 import haplotype_caller as haplotype_caller3
from lib.variants_freebayes import freebayes_caller, edit_freebayes_vcf
from lib.picard_actions import picard_sort
from lib.variants_octopus import octopus_caller
from lib.utils import move_bam, compile_identity
from lib.GATK_vcf import vcf_comparison
from lib.snpEff_ann import annotate_merged
from lib.picard_qc import get_coverage
from lib.picard_metrics import get_yield, get_hs_metrics, get_align_summary
from lib.extract_identity import mpileup, create_identity_table
from lib.process_identity import barcoding
from lib.count2 import extract_counts
from lib.variants_table import extract_info


# main configuration file
# couch_credentials = open('lib/config/couchdb').read().splitlines()
# run configuration file - should be moved to some argument parsing lib

# checks current version
VERSIONFILE = "VERSION"
VERSION = open(VERSIONFILE).read().strip()

logging.basicConfig(level=logging.INFO)


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
        if not Path(f"{datadir}/BAM/{sample}").exists():
            console.log(f"Creating {datadir}/BAM/{sample}")
            os.makedirs(f"{datadir}/BAM/{sample}")
        if not Path(f"{datadir}/BAM/{sample}/BAM/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/BAM/")
            os.makedirs(f"{datadir}/BAM/{sample}/BAM/")
        if not Path(f"{datadir}/BAM/{sample}/VCF/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/VCF/")
            os.makedirs(f"{datadir}/BAM/{sample}/VCF/")
        if not Path(f"{datadir}/BAM/{sample}/QC/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/QC/")
            os.makedirs(f"{datadir}/BAM/{sample}/QC/")
        if not Path(f"{datadir}/BAM/{sample}/Metrics/").exists():
            console.log(f"Creating {datadir}/BAM/{sample}/Metrics/")
            os.makedirs(f"{datadir}/BAM/{sample}/Metrics/")

    console.log("Directories created")


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

    for sample in samples:
        console.log(f"Processing {sample}")
        fastq_files = glob.glob(f"{datadir}/{sample}*.fastq.gz")
        console.log(f"Aligning sample {sample}")
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
    panel_samples = list(configuration["BED"].keys())
    fastqs = find_fastq(datadir, panel_samples, panel)
    sample_ids = get_ids(fastqs)

    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    console.log(f"Using BWA: {bwa}")
    console.log(f"Using SAMTOOLS: {samtools}")
    console.log(f"Using reference: {reference}")

    create_directories(datadir, sample_ids, panel)

    align_files(config, datadir, sample_ids, fastqs)

    analyse_pairs(config, datadir, sample_ids)

    return True


def analyse_pairs(config, datadir, samples):
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

        GATKvariants = haplotype_caller(
            datadir, sample, reference, bed_file[sample], gatk
        )

        GATK3variants = haplotype_caller3(
            datadir, sample, reference, bed_file[sample], gatk3
        )

        freebayesvariants = freebayes_caller(
            datadir, sample, reference, bed_file[sample], freebayes
        )
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

        if configuration["FinalDir"].find("Cplus") >= 1:
            to_return[sample]["cnv"] = extract_counts(
                datadir, "/apps/data/src/BED/new/C+_ALL_IDPE_01JUN2021_Window.bed", sample
            )
        else:
            to_return[sample]["cnv"] = extract_counts(
                datadir, "/apps/data/src/BED/new/CardiacALL_29MAR2021_Window.bed", sample
            )

        to_return[sample]["variants_table"] = extract_info(
            sample, datadir, configuration
        )
        console.log(f"Processing {sample} completed")



    if not Path(f"{datadir}/identity.txt").exists():
        console.log("Identity file does not exist, creating it")
        compile_identity(datadir)
    else:
        console.log("Identity file exists")



    #
    # # ########################################### #
    # #                                             #
    # #                  Identity                   #
    # #                                             #
    # # ########################################### #
    # try:
    #     if os.path.isfile(datadir + "/identity.txt"):
    #         logger.info("Identity file exists")
    #     else:
    #         logger.info("Compiling identity file")
    #         compile_identity(datadir)
    #         logger.info("Identity file compiled")
    # except Exception as e:
    #     logger.error("Some errors on identity determination, please check " + str(e))
    #
    # # ########################################### #
    # #                                             #
    # #                  Barcode                    #
    # #                                             #
    # # ########################################### #
    # logger.info("Barcode compilation")
    # try:
    #     process_identity.compile_barcodes(datadir)
    # except Exception as e:
    #     logger.error("Barcode compilation not possible" + str(e))
    #
    # # ########################################### #
    # #                                             #
    # #                  CNV calc                   #
    # #                                             #
    # # ########################################### #
    # logger.info("Calculating and saving CNV read normalization")
    # try:
    #     all_cnvs = cnv.compile_samples(datadir)
    #     cnv.cnv_calculation(datadir, all_cnvs, yaml_file)
    #     logger.info("Calculation completed, CNV mean file saved")
    # except Exception as e:
    #     logger.error("Some errors on CNV determination, please check " + str(e))
    #
    # # ########################################### #
    # #                                             #
    # #                Uniformity                   #
    # #                                             #
    # # ########################################### #
    # uniformity.get_coverage_values(datadir)
    #
    # # ########################################### #
    # #                                             #
    # #                Run DB creation              #
    # #                                             #
    # # ########################################### #
    #
    # process_folder(datadir)
    # add_samples.create_db_items(datadir, VERSION)
    #
    # # ########################################### #
    # #                                             #
    # #                  Variants DB                #
    # #                                             #
    # # ########################################### #
    # for pair in pairs:
    #     # try:
    #     #     to_return[pair]["final_report"] = generate_final_report.save_report(pair, datadir, transcript_location, config["Datadir"])
    #     #     try:
    #     #         copytree("resources", datadir + "/BAM/" + pair + "/resources")
    #     #     except Exception as e:
    #     #         logger.warning('Resources folder already exists ' + str(e))
    #     # except Exception as e:
    #     #     logger.error(str(e))
    #     # # check_sample.check_files(pair, datadir)
    #     check_sample.organize_files(pair, datadir)
    #
    # # generate main index deprecated
    # # generate_index.generate_index_page(datadir)
    #
    # # add metrics information to DB
    # try:
    #     # run_metrics_interop.insert_into_db(datadir)
    #     run_metrics_interop_couch.insert_into_db(datadir)
    # except Exception as e:
    #     logger.warning("Problems adding info to db " + str(e))
    #
    # # ########################################### #
    # #                                             #
    # #                  Samples DB                 #
    # #                                             #
    # # ########################################### #
    # try:
    #     add_samples.create_db_items(datadir)
    # except Exception as e:
    #     logger.warning("Problems adding samples to DB " + str(e))
    #
    # # ########################################### #
    # #                                             #
    # #                  Enrichment                 #
    # #                                             #
    # # ########################################### #
    # try:
    #     for pair in pairs:
    #         enrichment.get_enrichment(pair, datadir)
    # except Exception as e:
    #     logger.warning("Problems calculating enrichment" + str(e))
    #
    # # ########################################### #
    # #                                             #
    # #                  Gather Metrics             #
    # #                                             #
    # # ########################################### #
    # gather_metrics.process_all_metrics(datadir)

    return to_return


def compile_identity(datadir):
    """
    Function that reads samples' identity files
    and compiles them in a single file

    :param datadir: run location

    :type datadir: string

    :return: boolean
    """

    all_identity = open(datadir + "/identity.txt", "w")
    for filename in glob.glob(datadir + "/BAM/*/*"):
        if filename.find("identity.txt") >= 0:
            console.log(f"Found {filename}")
            all_identity.write(filename.split("/")[-2] + "\n")
            single_identity = open(filename).read()
            all_identity.write(single_identity)
    all_identity.close()

    return True


def generate_analysis(config, datadir, samples, panel):
    """
    Main function of the script, find input files, start alignemt and processing

    :param datadir: location of the FASTQ files

    :type datadir: string

    :return: patients dictionary
    :rtype: dictionary
    """

    console.log("Starting analysis and qc report generation")
    console.log(f"Configuration file: {config}")
    console.log(f"Datadir: {datadir}")
    if samples == []:
        console.log("Samples: ALL")
    else:
        console.log(f"Samples: {' '.join(samples)}")

    s = process_dir(config, datadir, samples, panel)

    # patients = ""
    # return patients


@click.command()
@click.option(
    "-c", "--configuration-file", "configuration_file", help="YAML file", required=True
)
@click.option(
    "-s",
    "--sample",
    "samples",
    help="sample to be analysed (it can be multiple, each with a -s)",
    multiple=True,
    required=False,
)
@click.option("-d", "--datadir", "datadir", help="run directory", required=True)
@click.option("-p", "--panel", "panel", help="panel to be used", required=True)
def run_analysis(configuration_file, samples, datadir, panel):
    """

    :param configuration_file:
    :param samples:
    :param datadir:
    :return:
    """

    samples = list(samples)

    console.log("Pipeline current version is " + VERSION)
    console.log("All requirements found, starting analysis")

    sample_dict = generate_analysis(configuration_file, datadir, samples, panel)


if __name__ == "__main__":
    run_analysis()

# 23GN-135G00063
# 23GN-135G00072
# 23GN-125G00029
