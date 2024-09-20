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
from typing import Dict, List, Optional
from datetime import datetime

import click
import yaml
from dotenv import dotenv_values
from rich.console import Console
from tinydb import TinyDB, Query

from lib import (
    bwa_align, cnv, count2, dup_indels, enrichment, extract_identity,
    GATK_vcf, picard_actions, picard_metrics, picard_qc, process_identity,
    recalibration, snpEff_ann, uniformity, utils, variants_freebayes,
    variants_GATK, variants_GATK3, variants_octopus
)

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
VERSIONFILE = Path(__file__).parent / "VERSION"
VERSION = VERSIONFILE.read_text().strip()

console = Console()


def get_db(datadir: Path) -> TinyDB:
    db_path = datadir / f"{datadir.name}_pipeline_logs.json"
    return TinyDB(db_path)


def log_to_db(db: TinyDB, message: str, level: str, program: str, sample_id: str, run_id: str):
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "message": message,
        "level": level,
        "program": program,
        "sample_id": sample_id,
        "run_id": run_id
    }
    db.insert(log_entry)


def split_string(ctx: object, param: object, value: Optional[str]) -> List[str]:
    if value is None:
        return []
    else:
        return [item.strip() for item in value.split(",")]


def find_fastq(datadir: Path, panel_samples: List[str], panel: str, db: TinyDB) -> List[Path]:
    fastqs = sorted(datadir.glob("*.fastq.gz"))
    for fastq in fastqs:
        message = f"Found {fastq}"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
    return fastqs


def get_ids(fastqs: List[Path]) -> List[str]:
    sample_ids = []
    for fastq in fastqs:
        sample_ids.append(Path(fastq).name.split("_")[0])
    return sorted(set(sample_ids))


def create_directories(datadir: Path, sample_ids: List[str], panel: str, db: TinyDB) -> None:
    try:
        (datadir / "BAM").mkdir(exist_ok=True)
    except FileExistsError:
        message = f"{datadir}/BAM/ already exists"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    message = f"Creating directories in {panel} directory"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    for sample in sample_ids:
        message = f"Processing {sample}"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", sample, datadir.name)

        for sub_dir in ["", "BAM", "VCF", "QC", "Metrics"]:
            path = datadir / "BAM" / sample / sub_dir
            path.mkdir(parents=True, exist_ok=True)
            message = f"Creating {path}"
            console.log(message)
            log_to_db(db, message, "INFO", "pipeline", sample, datadir.name)

    message = "Directories created/existed"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)


def align_files(config: Path, datadir: Path, samples: List[str], fastqs: List[Path], db: TinyDB) -> bool:
    env = dotenv_values(Path.cwd() / ".env")
    configuration = yaml.safe_load(config.read_text())
    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    message = "Aligning files"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    for sample in samples:
        message = f"Processing {sample}"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", sample, datadir.name)

        fastq_files = [fastq for fastq in fastqs if sample in fastq.name]
        message = f"Aligning sample {sample}"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", sample, datadir.name)

        run_bwa(sample, fastq_files, datadir, reference, bwa, samtools)

    return True


def process_dir(config: Path, datadir: Path, samples: List[str], panel: str, full_analysis: bool, db: TinyDB) -> Dict[
    str, Dict[str, bool]]:
    env = dotenv_values(Path.cwd() / ".env")
    configuration = yaml.safe_load(config.read_text())

    message = "Looking for FASTQ files"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    panel_samples = list(configuration["BED"].keys())

    fastqs = find_fastq(datadir, panel_samples, panel, db)
    sample_ids = get_ids(fastqs)

    message = f"Found {len(sample_ids)} samples"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    if len(samples) >= 1:
        message = f"Samples to be analysed: {' '.join(samples)}"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
        sample_ids = samples

    bwa = env["BWA"]
    samtools = env["SAMTOOLS"]
    reference = configuration["Reference"]

    message = f"Using BWA: {bwa}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    message = f"Using SAMTOOLS: {samtools}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    message = f"Using reference: {reference}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    create_directories(datadir, sample_ids, panel, db)
    align_files(config, datadir, sample_ids, fastqs, db)
    analyse_pairs(config, datadir, sample_ids, panel, full_analysis, db)

    return True


def analyse_pairs(config: Path, datadir: Path, samples: List[str], panel: str, full_analysis: bool, db: TinyDB) -> Dict[
    str, Dict[str, bool]]:
    to_return = defaultdict(dict)
    env = dotenv_values(Path.cwd() / ".env")
    configuration = yaml.safe_load(config.read_text())

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
        message = f"Processing {sample} :: {pos + 1} of {len(samples)}"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", sample, datadir.name)

        recalibration_step1 = base_recal1(datadir, sample, bed_file[sample], vcf_file, reference, gatk)
        to_return[sample]["recalibration1"] = recalibration_step1

        recalibration_final = recalibrate(datadir, sample, reference, gatk)
        to_return[sample]["recalibrate"] = recalibration_final
        move_bam(datadir, sample, "recal_reads")

        haplotype_caller(datadir, sample, reference, bed_file[sample], gatk)
        haplotype_caller3(datadir, sample, reference, bed_file[sample], gatk3)
        freebayes_caller(datadir, sample, reference, bed_file[sample], freebayes)
        to_return[sample]["picard_sort"] = picard_sort(datadir, sample, reference, picard)
        to_return[sample]["freebayes_edit"] = edit_freebayes_vcf(sample, datadir)
        to_return[sample]["variants_octopus"] = octopus_caller(datadir, sample, reference, bed_file[sample], octopus)
        to_return[sample]["vcf_merge"] = vcf_comparison(datadir, sample, reference, gatk3)
        to_return[sample]["snpEff"] = annotate_merged(sample, datadir, snpEff)
        to_return[sample]["picard_coverage"] = get_coverage(sample, datadir, reference, bait_file, picard)
        to_return[sample]["picard_coverage_panel"] = get_coverage(sample, datadir, reference, bed_file[sample], picard,
                                                                  "panel")
        to_return[sample]["picard_yield"] = get_yield(sample, datadir, picard)
        to_return[sample]["picard_hs_metrics"] = get_hs_metrics(sample, datadir, reference, bait_file, picard)
        to_return[sample]["picard_hs_metrics_panel"] = get_hs_metrics(sample, datadir, reference, bed_file[sample],
                                                                      picard, "panel")
        to_return[sample]["picard_align_metrics"] = get_align_summary(sample, datadir, reference, picard)
        to_return[sample]["mpileup_ident"] = mpileup(sample, datadir, "/apps/data/src/bundle/identity.txt", samtools)
        to_return[sample]["identity_table"] = create_identity_table(sample, datadir)
        to_return[sample]["full_identity"] = barcoding(sample, datadir)

        if panel == "Cplus":
            to_return[sample]["cnv"] = extract_counts(datadir,
                                                      "/apps/data/src/BED/new/C+_ALL_IDPE_01JUN2021_Window.bed", sample)
        else:
            to_return[sample]["cnv"] = extract_counts(datadir, "/apps/data/src/BED/new/CardiacALL_29MAR2021_Window.bed",
                                                      sample)

        message = f"Processing {sample} completed"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", sample, datadir.name)

    message = "Calculating uniformity"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
    get_coverage_values(datadir, panel)

    message = "Calculating enrichment"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
    for pos, sample in enumerate(samples):
        get_enrichment(sample, datadir, panel)

    if full_analysis:
        message = "Compiling identity file"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
        if not (datadir / "identity.txt").exists():
            message = "Identity file does not exist, creating it"
            console.log(message)
            log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
            compile_identity(datadir)
        else:
            message = "Identity file exists"
            console.log(message)
            log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

        message = "Compiling barcodes"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
        compile_barcodes(datadir)

        message = "Calculating and saving CNV read normalization"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
        all_cnvs = compile_samples(datadir)
        cnv_calculation(datadir, all_cnvs, config)
    else:
        message = "Full analysis not requested, skipping additional steps"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    return to_return


def generate_analysis(config: Path, datadir: Path, samples: List[str], panel: str, full_analysis: bool, db: TinyDB) -> Dict[str, Dict[str, bool]]:
    message = "Starting analysis and qc report generation"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    message = f"Configuration file: {config}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    message = f"Datadir: {datadir}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    message = f"Panel: {panel}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    if not samples:
        message = "No samples provided, will analyse all samples in the run"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)
        samples = []

    message = f"Samples: {' '.join(samples)}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir.name)

    s = process_dir(config, datadir, samples, panel, full_analysis, db)

    return s


@click.command
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
@click.option(
    "-f", "--full_analysis", "full_analysis", help="full analysis", is_flag=True
)
def run_analysis(
    configuration_file: str,
    datadir: str,
    panel: str,
    samples: List[str],
    full_analysis: bool,
) -> Dict[str, Dict[str, bool]]:

    config_path = Path(configuration_file)
    datadir_path = Path(datadir)
    db = get_db(datadir_path)

    if not samples:
        message = "No samples provided, will analyse all samples in the run"
        console.log(message)
        log_to_db(db, message, "INFO", "pipeline", "NA", datadir_path.name)
        samples = []

    message = f"Pipeline current version is {VERSION}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir_path.name)

    message = "All requirements found, starting analysis"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir_path.name)

    message = f"Full analysis: {full_analysis}"
    console.log(message)
    log_to_db(db, message, "INFO", "pipeline", "NA", datadir_path.name)

    sample_dict = generate_analysis(
        config_path, datadir_path, samples, panel, full_analysis, db
    )

    return sample_dict


if __name__ == "__main__":
    run_analysis()
