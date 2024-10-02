"""
.. module:: process_identity
    :platform: any
    :synopsis: Module that processes identity files and returns a full table with barcodes and other info
.. moduleauthor:: Paulo Nuin, March 2018, Updated by Assistant, October 2024
"""

import glob
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn
from rich.table import Table

from .log_api import log_to_api
from .db_logger import log_to_db
from tinydb import TinyDB

console = Console()


# positions for all genes used
GENES = {
    29416481: "ALK_E29",
    215645464: "BARD1_E4",
    14187449: "XPC_E16",
    112164561: "APC_E16",
    55214348: "EGFR_E4",
    116435768: "MET_E20",
    116436022: "MET_E21",
    30999280: "WRN_E26",
    90967711: "NBN_E10",
    145742879: "RECQL4_E3",
    98209594: "PTCH1_E23",
    133208979: "POLE_E45",
    32936646: "BRCA2_E17",
    89858602: "FANCI_E37",
    89849480: "FANCA_E16",
    59763347: "BRIP1_E19",
}

# hgvs notation for all variants/position
HGVS = {
    29416481: "c.4472A>G",
    215645464: "c.1134G>C",
    14187449: "c.2815C>A",
    112164561: "c.1635G>A",
    55214348: "c.474C>T",
    116435768: "c.3912C>T",
    116436022: "c.4071G>A",
    30999280: "c.3222G>T",
    90967711: "c.1197T>C",
    145742879: "c.132A>G",
    98209594: "c.3944C>T",
    133208979: "c.6252A>G",
    32936646: "c.7806-14T>C",
    89858602: "c.3906T>C",
    89849480: "c.1501G>A",
    59763347: "c.2755T>C",
}

# codes used in the identity report
CODES = {
    "ALK_E29000T": 1,
    "ALK_E290C0T": 2,
    "ALK_E290C00": 3,
    "ALK_E29": 5,
    "BARD1_E40C00": 1,
    "BARD1_E40CG0": 2,
    "BARD1_E400G0": 3,
    "BARD1_E4": 5,
    "XPC_E1600G0": 1,
    "XPC_E1600GT": 2,
    "XPC_E16000T": 3,
    "XPC_E16": 5,
    "APC_E1600G0": 1,
    "APC_E16A0G0": 2,
    "APC_E16A000": 3,
    "APC_E16": 5,
    "EGFR_E40C00": 1,
    "EGFR_E40C0T": 2,
    "EGFR_E4000T": 3,
    "EGFR_E4": 5,
    "MET_E200C00": 1,
    "MET_E200C0T": 2,
    "MET_E20000T": 3,
    "MET_E20": 5,
    "MET_E2100G0": 1,
    "MET_E21A0G0": 2,
    "MET_E21A000": 3,
    "MET_E21": 5,
    "WRN_E2600G0": 1,
    "WRN_E2600GT": 2,
    "WRN_E26000T": 3,
    "WRN_E26": 5,
    "NBN_E10A000": 1,
    "NBN_E10A0G0": 2,
    "NBN_E1000G0": 3,
    "NBN_E10": 5,
    "RECQL4_E3000T": 1,
    "RECQL4_E30C0T": 2,
    "RECQL4_E30C00": 3,
    "RECQL4_E3": 5,
    "PTCH1_E2300G0": 1,
    "PTCH1_E23A0G0": 2,
    "PTCH1_E23A000": 3,
    "PTCH1_E23": 5,
    "POLE_E45000T": 1,
    "POLE_E450C0T": 2,
    "POLE_E450C00": 3,
    "POLE_E45": 5,
    "BRCA2_E17000T": 1,
    "BRCA2_E170C0T": 2,
    "BRCA2_E170C00": 3,
    "BRCA2_E17": 5,
    "FANCI_E37000T": 1,
    "FANCI_E370C0T": 2,
    "FANCI_E370C00": 3,
    "FANCI_E37": 5,
    "FANCA_E160C00": 1,
    "FANCA_E160C0T": 2,
    "FANCA_E16000T": 3,
    "FANCA_E16": 5,
    "BRIP1_E19A000": 1,
    "BRIP1_E19A0G0": 2,
    "BRIP1_E1900G0": 3,
    "BRIP1_E19": 5,
    "MET_E210CG0": "U",
}





def read_identity(sample_id: str, datadir: Path, db: TinyDB) -> pd.DataFrame:
    """
    Read the full identity for a sample and return a DataFrame with information.

    :param sample_id: ID of the sample
    :param datadir: Run location
    :param db: Database for logging
    :return: DataFrame with sample data
    """
    identity_file = datadir / "BAM" / sample_id / "identity.txt"
    log_to_db(db, f"Reading identity file for {sample_id}", "INFO", "read_identity", sample_id, datadir.name)

    try:
        sample_identity = pd.read_csv(
            identity_file,
            sep="\t",
            names=["Chromosome", "Position", "A", "C", "G", "T"]
        )
        log_to_db(db, f"Successfully read identity file for {sample_id}", "INFO", "read_identity", sample_id,
                  datadir.name)
        return sample_identity
    except Exception as e:
        error_msg = f"Error reading identity file for {sample_id}: {str(e)}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_db(db, error_msg, "ERROR", "read_identity", sample_id, datadir.name)
        raise


def generate_barcode(sample_id: str, datadir: Path, sample_identity: pd.DataFrame, db: TinyDB) -> None:
    """
    Generate a numeric barcode for a sample and save it to a file.

    :param sample_id: ID of the sample
    :param datadir: Run location
    :param sample_identity: DataFrame with nucleotides in each position of the barcode
    :param db: Database for logging
    """
    barcode_order = [0, 3, 1, 12, 15, 4, 14, 13, 5, 6, 8, 11, 10, 9, 7, 2]
    barcode = ''.join(str(CODES[sample_identity["Full code"][i]]) for i in barcode_order)

    barcode_file = datadir / "BAM" / sample_id / f"{sample_id}.barcode.txt"

    try:
        barcode_file.write_text(barcode)
        console.print(Panel(f"[green]Barcode for {sample_id} is {barcode}[/green]"))
        log_to_db(db, f"Generated barcode for {sample_id}: {barcode}", "INFO", "generate_barcode", sample_id,
                  datadir.name)
    except Exception as e:
        error_msg = f"Error writing barcode file for {sample_id}: {str(e)}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_db(db, error_msg, "ERROR", "generate_barcode", sample_id, datadir.name)
        raise


def process_identity(sample_id: str, datadir: Path, sample_identity: pd.DataFrame, db: TinyDB) -> pd.DataFrame:
    """
    Process the identity for a sample and return a DataFrame with information.

    :param sample_id: ID of the sample
    :param datadir: Run location
    :param sample_identity: DataFrame with sample identity information
    :param db: Database for logging
    :return: DataFrame with processed sample data
    """
    log_to_db(db, f"Processing identity for {sample_id}", "INFO", "process_identity", sample_id, datadir.name)

    try:
        nucleotides = ["A", "C", "G", "T"]

        sample_identity["ID"] = sample_identity["Position"].map(GENES)
        sample_identity["HGVS"] = sample_identity["Position"].map(HGVS)
        sample_identity["Total reads"] = sample_identity[nucleotides].sum(axis=1)

        for nucleotide in nucleotides:
            sample_identity[f"pc{nucleotide}"] = sample_identity[nucleotide] / sample_identity["Total reads"]
            sample_identity[f"code{nucleotide}"] = np.where(sample_identity[f"pc{nucleotide}"] > 0.1, nucleotide, "0")

        sample_identity["Full code"] = sample_identity["ID"] + sample_identity["codeA"] + sample_identity["codeC"] + \
                                       sample_identity["codeG"] + sample_identity["codeT"]

        output_file = datadir / "BAM" / sample_id / f"{sample_id}.identity_full.txt"
        sample_identity.to_csv(output_file, sep="\t", index=False)

        console.print(Panel(f"[green]Full identity file generated for {sample_id}[/green]"))
        log_to_db(db, f"Full identity file generated for {sample_id}", "INFO", "process_identity", sample_id,
                  datadir.name)

        return sample_identity
    except Exception as e:
        error_msg = f"Error processing identity for {sample_id}: {str(e)}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_db(db, error_msg, "ERROR", "process_identity", sample_id, datadir.name)
        raise


def barcoding(sample_id: str, datadir: Path, db: TinyDB) -> None:
    """
    Process the identity for a sample, generate a numeric barcode, and save the results.

    :param sample_id: ID of the sample
    :param datadir: Run location
    :param db: Database for logging
    """
    console.print(Panel(f"[bold blue]Processing identity and generating barcode for {sample_id}[/bold blue]"))
    log_to_db(db, f"Starting identity processing and barcode generation for {sample_id}", "INFO", "barcoding",
              sample_id, datadir.name)

    try:
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn()) as progress:
            task1 = progress.add_task("[green]Reading identity...", total=100)
            sample_identity = read_identity(sample_id, datadir, db)
            progress.update(task1, advance=100)

            task2 = progress.add_task("[cyan]Processing identity...", total=100)
            processed_identity = process_identity(sample_id, datadir, sample_identity, db)
            progress.update(task2, advance=100)

            task3 = progress.add_task("[magenta]Generating barcode...", total=100)
            generate_barcode(sample_id, datadir, processed_identity, db)
            progress.update(task3, advance=100)

        console.print(
            Panel(f"[bold green]Successfully processed identity and generated barcode for {sample_id}[/bold green]"))
        log_to_db(db, f"Successfully processed identity and generated barcode for {sample_id}", "INFO", "barcoding",
                  sample_id, datadir.name)
    except Exception as e:
        console.print(
            Panel(f"[bold red]Error processing identity and generating barcode for {sample_id}: {str(e)}[/bold red]"))
        log_to_db(db, f"Error processing identity and generating barcode for {sample_id}: {str(e)}", "ERROR",
                  "barcoding", sample_id, datadir.name)
        raise


def find_duplicates(datadir: Path, db: TinyDB) -> Dict[str, bool]:
    """
    Check for duplicate barcodes in the data directory.

    :param datadir: The directory where the data is located
    :param db: Database for logging
    :return: A dictionary where the keys are sample IDs and the values are boolean indicating whether the barcode is duplicated
    """
    log_to_db(db, "Checking for duplicate barcodes", "INFO", "find_duplicates", "NA", datadir.name)

    try:
        barcodes_file = datadir / "barcodes.txt"
        with barcodes_file.open() as f:
            barcodes = f.read().splitlines()

        samples_dict = {line.split("\t")[0]: line.split("\t")[1] for line in barcodes}
        codes = list(samples_dict.values())

        check = {sample: codes.count(code) > 1 for sample, code in samples_dict.items()}

        log_to_db(db, "Duplicate barcode check completed", "INFO", "find_duplicates", "NA", datadir.name)
        return check
    except Exception as e:
        error_msg = f"Error checking for duplicate barcodes: {str(e)}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_db(db, error_msg, "ERROR", "find_duplicates", "NA", datadir.name)
        raise


def compile_barcodes(datadir: Path, db: TinyDB) -> str:
    """
    Compile all barcodes from the data directory into a single file.

    :param datadir: The directory where the data is located
    :param db: Database for logging
    :return: A string indicating the success of the operation
    """
    console.print(Panel("[bold blue]Compiling barcodes[/bold blue]"))
    log_to_db(db, "Starting barcode compilation", "INFO", "compile_barcodes", "NA", datadir.name)

    try:
        barcodes_file = datadir / "barcodes.txt"
        barcode_dict = {}

        for sample_dir in (datadir / "BAM").glob("*"):
            sample_id = sample_dir.name
            barcode_file = sample_dir / f"{sample_id}.barcode.txt"

            if barcode_file.exists():
                barcode = barcode_file.read_text().strip()
                barcode_dict[sample_id] = barcode

        with barcodes_file.open("w") as f:
            for sample_id, barcode in barcode_dict.items():
                f.write(f"{sample_id}\t{barcode}\n")

        console.print(Panel("[bold green]Barcodes compiled successfully[/bold green]"))
        log_to_db(db, "Barcodes compiled successfully", "INFO", "compile_barcodes", "NA", datadir.name)
        return "success"
    except Exception as e:
        error_msg = f"Error compiling barcodes: {str(e)}"
        console.print(Panel(f"[bold red]{error_msg}[/bold red]"))
        log_to_db(db, error_msg, "ERROR", "compile_barcodes", "NA", datadir.name)
        raise


def process_barcodes(datadir: Path, samples: List[str], db: TinyDB) -> None:
    """
    Process barcodes for all samples in the data directory.

    :param datadir: The directory where the data is located
    :param samples: List of sample IDs to process
    :param db: Database for logging
    """
    console.print(Panel("[bold blue]Processing barcodes for all samples[/bold blue]"))
    log_to_db(db, "Starting barcode processing for all samples", "INFO", "process_barcodes", "NA", datadir.name)

    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn()) as progress:
        task = progress.add_task("[cyan]Processing samples...", total=len(samples))

        for sample in samples:
            try:
                barcoding(sample, datadir, db)
                progress.update(task, advance=1)
            except Exception as e:
                console.print(f"[bold red]Error processing barcode for {sample}: {str(e)}[/bold red]")
                log_to_db(db, f"Error processing barcode for {sample}: {str(e)}", "ERROR", "process_barcodes", sample,
                          datadir.name)

    compile_barcodes(datadir, db)
    duplicate_check = find_duplicates(datadir, db)

    # Display duplicate check results
    table = Table(title="Duplicate Barcode Check")
    table.add_column("Sample ID", style="cyan")
    table.add_column("Is Duplicate", style="magenta")

    for sample, is_duplicate in duplicate_check.items():
        table.add_row(sample, str(is_duplicate))

    console.print(table)

    console.print(Panel("[bold green]Finished processing barcodes for all samples[/bold green]"))
    log_to_db(db, "Finished barcode processing for all samples", "INFO", "process_barcodes", "NA", datadir.name)


if __name__ == "__main__":
    # Example usage
    from tinydb import TinyDB

    datadir = Path("/path/to/data")
    db = TinyDB(datadir / "pipeline_logs.json")
    samples = ["sample1", "sample2", "sample3"]

    process_barcodes(datadir, samples, db)