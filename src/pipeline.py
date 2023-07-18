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
from pathlib import Path, PurePath
import click
import yaml

from rich.console import Console

console = Console()

from dotenv import dotenv_values

from lib.bwa_align import run_bwa

# main configuration file
# couch_credentials = open('lib/config/couchdb').read().splitlines()
# run configuration file - should be moved to some argument parsing lib

# checks current version
VERSIONFILE = "VERSION"
VERSION = open(VERSIONFILE).read().strip()

logging.basicConfig(level=logging.INFO)


def get_code(sample_id):
    """
    Function that returns the  BED code of the sample
    :param sample_id:
    :return:
    """

    return sample_id[-2:]


def process_dir(config, datadir, samples):
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

    # call the function to find the FASTQ files

    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    console.log(f"Using BWA: {bwa}")
    console.log(f"Using SAMTOOLS: {samtools}")
    console.log(f"Using reference: {reference}")

    # return variable - should be used in logging or some centralized action. queueing
    to_return = defaultdict(dict)

    # align all pairs
    for sample in samples:
        console.log(f"Processing {sample}")
        patient_bwa = {}
        if not os.path.exists(f"{datadir}/BAM/{sample}/BAM/"):
            console.log(f"Creating {datadir}/BAM/{sample}")
            os.makedirs(f"{datadir}/BAM/{sample}/BAM/")
        if not os.path.exists(f"{datadir}/BAM/{sample}/VCF/"):
            console.log(f"Creating {datadir}/BAM/{sample}/VCF/")
            os.makedirs(f"{datadir}/BAM/{sample}/VCF/")
        if not os.path.exists(f"{datadir}/BAM/{sample}/QC/"):
            console.log(f"Creating {datadir}/BAM/{sample}/QC/")
            os.makedirs(f"{datadir}/BAM/{sample}/QC/")
        if not os.path.exists(f"{datadir}/BAM/{sample}/Metrics/"):
            console.log(f"Creating {datadir}/BAM/{sample}/Metrics/")
            os.makedirs(f"{datadir}/BAM/{sample}/Metrics/")

        # code = get_code(pair)
        # create directories in sample locations
        # alignment
        fastqs = glob.glob(f"{datadir}/Data/Intensities/BaseCalls/{sample}*.fastq.gz")
        console.log(f"Found {len(fastqs)} FASTQ files\n {' '.join(fastqs)}")
        to_return[sample]["bam"] = run_bwa(
            sample, fastqs, datadir, reference, bwa, samtools
        )
        # patient_bwa[pair] = ("BWA complete", str(datetime.now()))
    #
    to_return.update(analyse_pairs(configuration, datadir, samples))
    #
    # return to_return


def parse_picard(data_lines, pair, picard_metric):
    """
    Function that parses picard different metrics files
    and updates the database with results

    :param data_lines: content of the picard metric file
    :param pair: sample ID
    :param picard_metric: type of metric contained in the file

    :type data_lines: string
    :type pair: string
    :type picard_metric: string

    :return: success
    :rtype: string
    """

    return "success"


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

    print(env)
    # load software configuration

    samtools = env["SAMTOOLS"]
    picard = env["PICARD"]
    gatk = env["GATK"]
    gatk3 = env["GATK3"]
    freebayes = env["FREEBAYES"]
    qualimap = env["QUALIMAP"]
    snpEff = env["SNPEFF"]
    annovar_dir = env["ANNOVAR"]
    # transcript_location = software_conf["transcripts"]
    octopus = env["OCTOPUS"]
    reference = configuration["Reference"]
    vcf_file = configuration["VCF"]
    bed_file = configuration["BED"]
    bait_file = configuration["BAIT"]


    # index = 0  # keeps track of the current sample
    # sorted_pairs = sorted(pairs.keys())
    #
    # # main analysis loop
    # for pair in sorted_pairs:
    #     try:
    #         # checks if sample is fully analysed before starting
    #         sample_status = check_sample.check_files(pair, datadir)
    #         index += 1
    #         logger.info(sample_status)
    #         sample_status = "Proceed"
    #         if sample_status == "Proceed":
    #             logger.info(
    #                 "Starting analysis of sample "
    #                 + pair
    #                 + " "
    #                 + str(index)
    #                 + " of "
    #                 + str(len(pairs))
    #             )
    #
    #             # remove_duplicates = dup_indels.samtools_duplicates(pair, datadir, samtools)
    #             remove_duplicates = dup_indels.remove_duplicates(pair, datadir, picard)
    #             to_return[pair]["dedup"] = remove_duplicates
    #
    #             # add groups to the BAM file to avoid errors
    #             adding_groups = dup_indels.add_groups(pair, datadir, picard, samtools)
    #             to_return[pair]["add_groups"] = adding_groups
    #             update_dict[pair] = ("Realign complete", str(datetime.now()))
    #
    #             # ########################################### #
    #             #                                             #
    #             #           Recalibration                     #
    #             #                                             #
    #             # ########################################### #
    #             # Recalibration step 1
    #             recalibration_step1 = recalibration.base_recal1(
    #                 pair, datadir, bed_file[pair], vcf_file, reference, gatk
    #             )
    #             to_return[pair]["recalibration1"] = recalibration_step1
    #
    #             # Recalibration step 2
    #             recalibration_final = recalibration.recalibrate(
    #                 pair, datadir, reference, gatk
    #             )
    #             to_return[pair]["recalibrate"] = recalibration_final
    #             update_dict2[pair] = ("Recalibration complete", str(datetime.now()))
    #
    #             # ########################################### #
    #             #                                             #
    #             #           Variant Calling                   #
    #             #                                             #
    #             # ########################################### #
    #             # GATK4
    #             GATKvariants = variants_GATK.haplotype_caller(
    #                 pair, datadir, reference, bed_file[pair], gatk
    #             )
    #
    #             # GATK3
    #             GATK3variants = variants_GATK3.haplotype_caller(
    #                 pair, datadir, reference, bed_file[pair], gatk3
    #             )
    #             to_return[pair]["variants_GATK_1"] = GATKvariants
    #             update_dict3[pair] = ("GATK VCF generated", str(datetime.now()))
    #
    #             # Freebayes
    #             freebayesvariants = variants_freebayes.freebayes(
    #                 pair, datadir, reference, bed_file[pair], freebayes
    #             )
    #             if freebayesvariants == "error":
    #                 to_return[pair]["variants_freebayes"] = freebayesvariants
    #                 raise ValueError("\t\tFreebayes variant called failed")
    #             else:
    #                 to_return[pair]["variants_freebayes"] = freebayesvariants
    #             update_dict4[pair] = ("Freebayes VCF generated", str(datetime.now()))
    #
    #             # Freebayes VCF modifications
    #             to_return[pair]["picard_sort"] = picard_actions.picard_sort(
    #                 pair, datadir, reference, picard
    #             )
    #             to_return[pair][
    #                 "freebayes_edit"
    #             ] = variants_freebayes.edit_freebayes_vcf(pair, datadir)
    #
    #             # Octopus
    #             to_return[pair]["variants_octopus"] = variants_octopus.octopus(
    #                 pair, datadir, reference, bed_file[pair], octopus
    #             )
    #
    #             # Merging VCFs
    #             to_return[pair]["vcf_merge"] = GATK_vcf.vcf_comparison(
    #                 pair, datadir, reference, gatk3
    #             )
    #
    #             # ########################################### #
    #             #                                             #
    #             #               Annotation                    #
    #             #                                             #
    #             # ########################################### #
    #             # snpEff
    #             to_return[pair]["snpEff"] = snpEff_ann.annotate_merged(
    #                 pair, datadir, snpEff
    #             )
    #             update_dict5[pair] = (
    #                 "snpEff annotation generated",
    #                 str(datetime.now()),
    #             )
    #
    #             to_return[pair]["snpEff_pseudo"] = snpEff_ann.annotate_pseudo(
    #                 pair, datadir, snpEff
    #             )
    #             update_dict5[pair] = (
    #                 "snpEff pseudo annotation generated",
    #                 str(datetime.now()),
    #             )
    #
    #             # Annovar
    #             to_return[pair][annovar] = annovar.get_annotation(
    #                 pair, datadir, annovar_dir
    #             )
    #             to_return[pair]["vcf_parser"] = vcf_parser.parse_vcf(
    #                 pair, datadir, transcript_location, config["FinalDir"]
    #             )
    #
    #             # ########################################### #
    #             #                                             #
    #             #                  QA/QC                      #
    #             #                                             #
    #             # ########################################### #
    #             # Qualimap
    #             to_return[pair]["qualimap"] = run_qualimap.call_qualimap(
    #                 pair, datadir, bed_file[pair], qualimap
    #             )
    #             update_dict6[pair].append("qualimap")
    #
    #             # Picard coverage
    #             to_return[pair]["picard_coverage"] = picard_qc.get_coverage(
    #                 pair, datadir, reference, bait_file, picard, transcript_location
    #             )
    #             to_return[pair]["picard_coverage_panel"] = picard_qc.get_coverage(
    #                 pair,
    #                 datadir,
    #                 reference,
    #                 bed_file[pair],
    #                 picard,
    #                 transcript_location,
    #                 "panel",
    #             )
    #
    #             update_dict6[pair].append("picard coverage")
    #
    #             to_return[pair]["picard_yield"] = picard_metrics.get_yield(
    #                 pair, datadir, picard
    #             )
    #             if os.path.isfile(datadir + "/BAM/" + pair + "/" + pair + ".yield.out"):
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir + "/BAM/" + pair + "/" + pair + ".yield.out"
    #                     ),
    #                     pair,
    #                     "yield",
    #                 )
    #             else:
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir + "/BAM/" + pair + "/Metrics/" + pair + ".yield.out"
    #                     ),
    #                     pair,
    #                     "yield",
    #                 )
    #             update_dict6[pair].append("picard yield")
    #
    #             to_return[pair]["picard_hs_metrics"] = picard_metrics.get_hs_metrics(
    #                 pair, datadir, reference, bait_file, picard
    #             )
    #             to_return[pair][
    #                 "picard_hs_metrics_panel"
    #             ] = picard_metrics.get_hs_metrics(
    #                 pair, datadir, reference, bed_file[pair], picard, "panel"
    #             )
    #
    #             if not os.path.isfile(
    #                 datadir + "/BAM/" + pair + "/Metrics/" + pair + ".nucl.out"
    #             ):
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir + "/BAM/" + pair + "/" + pair + ".hs_metrics.out"
    #                     ),
    #                     pair,
    #                     "hs",
    #                 )
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir
    #                         + "/BAM/"
    #                         + pair
    #                         + "/"
    #                         + pair
    #                         + ".hs_metrics.panel.out"
    #                     ),
    #                     pair,
    #                     "hs_panel",
    #                 )
    #             else:
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir
    #                         + "/BAM/"
    #                         + pair
    #                         + "/Metrics/"
    #                         + pair
    #                         + ".hs_metrics.out"
    #                     ),
    #                     pair,
    #                     "hs",
    #                 )
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir
    #                         + "/BAM/"
    #                         + pair
    #                         + "/Metrics/"
    #                         + pair
    #                         + ".hs_metrics.panel.out"
    #                     ),
    #                     pair,
    #                     "hs_panel",
    #                 )
    #
    #             update_dict6[pair].append("picard hs metrics")
    #             to_return[pair][
    #                 "picard_align_metrics"
    #             ] = picard_metrics.get_align_summary(pair, datadir, reference, picard)
    #
    #             if not os.path.isfile(
    #                 datadir + "/BAM/" + pair + "/Metrics/" + pair + ".nucl.out"
    #             ):
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir + "/BAM/" + pair + "/" + pair + ".align_metrics.out"
    #                     ),
    #                     pair,
    #                     "align",
    #                 )
    #             else:
    #                 parse_picard(
    #                     picard_parser.main_parser(
    #                         datadir
    #                         + "/BAM/"
    #                         + pair
    #                         + "/Metrics/"
    #                         + pair
    #                         + ".align_metrics.out"
    #                     ),
    #                     pair,
    #                     "align",
    #                 )
    #             update_dict6[pair].append("picard align metrics")
    #
    #             # ########################################### #
    #             #                                             #
    #             #                Identity.                    #
    #             #                                             #
    #             # ########################################### #
    #             to_return[pair]["mpileup_ident"] = extract_identity.mpileup(
    #                 pair, datadir, "/opt/bundle/identity.txt", samtools
    #             )
    #             to_return[pair][
    #                 "identity_table"
    #             ] = extract_identity.create_identity_table(pair, datadir)
    #             to_return[pair]["full_identity"] = process_identity.barcoding(
    #                 pair, datadir
    #             )
    #
    #             # ########################################### #
    #             #                                             #
    #             #                  Sample CNV                 #
    #             #                                             #
    #             # ########################################### #
    #             if config["FinalDir"].find("Cplus") >= 1:
    #                 to_return[pair]["cnv"] = count2.extract_counts(
    #                     datadir, "/opt/BED/new/C+_ALL_IDPE_01JUN2021_Window.bed", pair
    #                 )
    #             else:
    #                 to_return[pair]["cnv"] = count2.extract_counts(
    #                     datadir, "/opt/BED/new/CardiacALL_29MAR2021_Window.bed", pair
    #                 )
    #
    #             to_return[pair]["variants_table"] = variants_table.extract_info(
    #                 pair, datadir, transcript_location, config["Datadir"]
    #             )
    #             logger.info("Process of sample " + pair + " completed")
    #         else:
    #             logger.info("Process of sample %s completed previously" % (pair))
    #     except Exception as e:
    #         logger.error("Process of sample %s failed with error %s" % (pair, e))
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
            logger.info(filename)
            all_identity.write(filename.split("/")[-2] + "\n")
            single_identity = open(filename).read()
            all_identity.write(single_identity)
    all_identity.close()

    return True


def generate_analysis(config, datadir, samples):
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
    console.log(f"Samples: {' '.join(samples)}")

    patients = process_dir(config, datadir, samples)

    patients = ""
    return patients


def check_datadir(datadir, configuration_file):
    """
    Function that checled the run directory for FASTQ files

    :param datadir: run location
    :param configuration_file: run yaml file configuration

    :type datadir: string
    :type configuration_file: string
    """

    temp = datadir.split("/")

    if len(temp) > 1:
        return datadir

    if datadir in configuration_file:
        return os.path.dirname(configuration_file)

    return "FASTQ source not found, please check the configuration and restart the analysis"


@click.command()
@click.option("-c", "--configuration-file", "configuration_file", help="YAML file")
@click.option(
    "-s",
    "--sample",
    "samples",
    help="sample to be analysed (it can be multiple, each with a -s)",
    multiple=True,
)
@click.option("-d", "--datadir", "datadir", help="run directory")
def run_analysis(configuration_file, samples, datadir):
    """

    :param configuration_file:
    :param samples:
    :param datadir:
    :return:
    """

    samples = list(samples)

    console.log("Pipeline current version is " + VERSION)
    console.log("All requirements found, starting analysis")

    patient_dict = generate_analysis(configuration_file, datadir, samples)


if __name__ == "__main__":
    run_analysis()

# 23GN-135G00063
# 23GN-135G00072
# 23GN-125G00029