"""
.. module:: align
    :platform: any
    :synopsis: This module is the main piece in the NGS pipeline, performing
                FASTQ alignments and calling all other steps in the process
.. moduleauthor:: Paulo Nuin, December 2015

"""

import os
from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Optional
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


# checks current version
VERSIONFILE = "VERSION"
VERSION = open(VERSIONFILE).read().strip()

console = Console()


def split_string(ctx: object, param: object, value: Optional[str]) -> List[str]:
    """
    Splits a string into a list of items. If the input string is None, returns an empty list.
    Each item in the list is stripped of leading and trailing whitespace.

    Parameters:
    ctx: The context in which the function is called. Not used in this function.
    param: The parameter being processed. Not used in this function.
    value (str): The string to be split. If None, an empty list is returned.

    Returns:
    list: A list of items obtained by splitting the input string by comma.
    """

    # If value is None, return an empty list
    if value is None:
        return []
    else:
        # Split the value by comma and return a list of stripped items
        return [item.strip() for item in value.split(",")]


def find_fastq(datadir: str, panel_samples: List[str], panel: str) -> List[str]:
    """
    Finds all the FASTQ files in the given directory.

    Parameters:
    datadir (str): The directory in which to search for FASTQ files.
    panel_samples (list): The list of panel samples. Not used in this function.
    panel (str): The panel to be used. Not used in this function.

    Returns:
    list: A list of found FASTQ files.
    """

    # Use glob to find all files that match the pattern "*.fastq.gz"
    # fastqs = glob.glob(f"{datadir}/*.fastq.gz")
    fastqs = list(Path(datadir).glob("*.fastq.gz"))

    # Sort the found files
    fastqs = sorted(fastqs)

    to_return = []
    for fastq in fastqs:
        # Log each found file
        console.log(f"Found {fastq}")
        log_to_api(f"Found {fastq}", "INFO", "pipeline", "NA", Path(datadir).name)
        to_return.append(Path(fastq).name)

    # Return the list of found FASTQ files
    return fastqs


def get_ids(fastqs: List[Path]) -> List[str]:
    """
    Extracts the sample IDs from the list of FASTQ files.

    Parameters:
    fastqs (List[Path]): The list of FASTQ files.

    Returns:
    List[str]: A sorted list of unique sample IDs.
    """

    sample_ids = []
    for fastq in fastqs:
        # Split the file name on underscore and take the first part as the sample ID
        sample_ids.append(Path(fastq).name.split("_")[0])

    # Remove duplicates and sort the sample IDs
    sample_ids = sorted(set(sample_ids))

    return sample_ids


def create_directories(datadir: str, sample_ids: List[str], panel: str) -> None:
    """
    Creates necessary directories for each sample in the given directory.

    Parameters:
    datadir (str): The directory in which to create sample directories.
    sample_ids (List[str]): The list of sample IDs for which directories are to be created.
    panel (str): The panel to be used. Not used in this function.

    Returns:
    None
    """

    # Try to create a directory for BAM files
    try:
        os.makedirs(f"{datadir}/BAM/")
    except FileExistsError:
        # If the directory already exists, log a message and continue
        console.log(f"{datadir}/BAM/ already exists")

    console.log(f"Creating directories in {panel} directory")

    # For each sample ID, create necessary directories
    for sample in sample_ids:
        console.log(f"Processing {sample}")
        log_to_api(
            f"Processing {sample}", "INFO", "pipeline", sample, Path(datadir).name
        )

        # Create directories for each sample if they do not exist
        for sub_dir in ["", "BAM", "VCF", "QC", "Metrics"]:
            path = Path(datadir) / "BAM" / sample / sub_dir
            path.mkdir(parents=True, exist_ok=True)
            console.log(f"Creating {path}")
            log_to_api(
                f"Creating {path}",
                "INFO",
                "pipeline",
                sample,
                Path(datadir).name,
            )

    console.log("Directories created/existed")


def align_files(
    config: str, datadir: Path, samples: List[str], fastqs: List[Path]
) -> bool:
    """
    Aligns the FASTQ files for each sample using the BWA aligner.

    Parameters:
    config (str): The configuration file for the run.
    datadir (Path): The directory containing the FASTQ files.
    samples (List[str]): The list of samples to be analysed.
    fastqs (List[Path]): The list of FASTQ files.

    Returns:
    bool: True if the alignment was successful, False otherwise.
    """

    # Load the environment variables and configuration
    env = dotenv_values(f"{Path.cwd()}/.env")
    configuration = yaml.safe_load(open(config))
    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    # Log the start of the alignment process
    console.log("Aligning files")
    log_to_api("Aligning files", "INFO", "pipeline", "NA", Path(datadir).name)

    # For each sample, align the corresponding FASTQ files
    for sample in samples:
        console.log(f"Processing {sample}")
        log_to_api(
            f"Processing {sample}", "INFO", "pipeline", sample, Path(datadir).name
        )
        fastq_files = [fastq for fastq in fastqs if sample in fastq.name]
        console.log(f"Aligning sample {sample}")
        log_to_api(
            f"Aligning sample {sample}", "INFO", "pipeline", sample, Path(datadir).name
        )
        run_bwa(sample, fastq_files, datadir, reference, bwa, samtools)

    return True


def process_dir(
    config: str, datadir: str, samples: List[str], panel: str
) -> Dict[str, Dict[str, bool]]:
    """
    Main function that orchestrates the entire pipeline of processing the data.

    Parameters:
    config (str): The configuration file for the run.
    datadir (str): The directory with the location of the FASTQ files.
    samples (List[str]): The list of samples in the run.
    panel (str): The panel to be used.

    Returns:
    Dict[str, Dict[str, bool]]: A dictionary with the success/failure of each step for each sample.
    """

    # Load the environment variables and configuration
    env = dotenv_values(f"{Path.cwd()}/.env")
    configuration = yaml.safe_load(open(config))

    # Log the start of the FASTQ file search
    console.log("Looking for FASTQ files")
    log_to_api("Looking for FASTQ files", "INFO", "pipeline", "NA", Path(datadir).name)

    # Get the list of panel samples from the configuration
    panel_samples = list(configuration["BED"].keys())

    # Find the FASTQ files and get the sample IDs
    fastqs = find_fastq(datadir, panel_samples, panel)
    sample_ids = get_ids(fastqs)

    # Log the number of found samples
    console.log(f"Found {len(sample_ids)} samples")
    log_to_api(
        f"Found {len(sample_ids)} samples", "INFO", "pipeline", "NA", Path(datadir).name
    )

    # If specific samples are provided, use them instead
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

    # Load the software configuration
    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    # Log the software configuration
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

    # Create necessary directories, align files, and analyse pairs
    create_directories(datadir, sample_ids, panel)
    align_files(config, datadir, sample_ids, fastqs)
    analyse_pairs(config, datadir, sample_ids, panel)

    return True


def analyse_pairs(
    config: str, datadir: str, samples: List[str], panel: str
) -> Dict[str, Dict[str, bool]]:
    """
    Orchestrates the entire pipeline of processing the data for each sample.

    Parameters:
    config (str): The configuration file for the run.
    datadir (str): The directory with the location of the FASTQ files.
    samples (List[str]): The list of samples in the run.
    panel (str): The panel to be used.

    Returns:
    Dict[str, Dict[str, bool]]: A dictionary with the success/failure of each step for each sample.

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


def generate_analysis(
    config: str, datadir: str, samples: List[str], panel: str
) -> Dict[str, Dict[str, bool]]:
    """
    Main function of the script, find input files, start alignment and processing.

    Parameters:
    config (str): The configuration file for the run.
    datadir (str): The directory with the location of the FASTQ files.
    samples (List[str]): The list of samples in the run.
    panel (str): The panel to be used.

    Returns:
    Dict[str, Dict[str, bool]]: A dictionary with the success/failure of each step for each sample.
    """

    # Log the start of the analysis and quality control report generation process
    console.log("Starting analysis and qc report generation")
    log_to_api(
        "Starting analysis and qc report generation",
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )

    # Log the configuration file, data directory, and panel
    console.log(f"Configuration file: {config}")
    log_to_api(
        f"Configuration file: {config}", "INFO", "pipeline", "NA", Path(datadir).name
    )
    console.log(f"Datadir: {datadir}")
    log_to_api(f"Datadir: {datadir}", "INFO", "pipeline", "NA", Path(datadir).name)
    console.log(f"Panel: {panel}")
    log_to_api(f"Panel: {panel}", "INFO", "pipeline", "NA", Path(datadir).name)

    # If no specific samples are provided, log a message and set samples to an empty list
    if not samples:
        console.log("No samples provided, will analyse all samples in the run")
        samples = []

    # Log the samples to be analysed
    console.log(f"Samples: {' '.join(samples)}")
    log_to_api(
        f"Samples: {' '.join(samples)}",
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )

    # Start the analysis process by calling the process_dir function
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
def run_analysis(
    configuration_file: str, datadir: str, panel: str, samples: List[str]
) -> Dict[str, Dict[str, bool]]:
    """
    Starts the analysis process by calling the generate_analysis function.

    Parameters:
    configuration_file (str): The YAML configuration file for the run.
    datadir (str): The directory of the run.
    panel (str): The panel to be used.
    samples (List[str]): The list of samples to be analysed.

    Returns:
    Dict[str, Dict[str, bool]]: A dictionary with the success/failure of each step for each sample.
    """

    # If no specific samples are provided, log a message and set samples to an empty list
    if not samples:
        console.log("No samples provided, will analyse all samples in the run")
        samples = []

    # Log the current version of the pipeline
    console.log("Pipeline current version is " + VERSION)
    log_to_api(
        "Pipeline current version is " + VERSION,
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )

    # Log the start of the analysis
    console.log("All requirements found, starting analysis")
    log_to_api(
        "All requirements found, starting analysis",
        "INFO",
        "pipeline",
        "NA",
        Path(datadir).name,
    )

    # Start the analysis process by calling the generate_analysis function
    sample_dict = generate_analysis(configuration_file, datadir, samples, panel)

    return sample_dict


if __name__ == "__main__":
    run_analysis()
