"""
.. module:: count2.py
    :platform: any
    :synopsis: Module to generate read counts under each exon
.. moduleauthor:: Paulo Nuin, October 2017, Updated September 2024
"""

from pathlib import Path
import pysam
from rich.console import Console
from rich.progress import Progress
from typing import Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
from .db_logger import log_to_db
from clickhouse_driver import Client
from dotenv import load_dotenv
import os
from datetime import datetime

import click
from pathlib import Path
from .db_logger import get_sample_db


console = Console()

# Load environment variables
load_dotenv()

# Access environment variables
SAMTOOLS = os.getenv('SAMTOOLS')
BED_FILES = os.getenv('BED_FILES')
DB_HOST = os.getenv('DB_HOST')
DB_PORT = int(os.getenv('DB_PORT', 9000))
DB_USER = os.getenv('DB_USER')
DB_PASSWORD = os.getenv('DB_PASSWORD')
DB_NAME = os.getenv('DB_NAME')

def get_clickhouse_client():
    return Client(
        host=DB_HOST,
        port=DB_PORT,
        user=DB_USER,
        password=DB_PASSWORD,
        database=DB_NAME
    )

def process_bed_region(samfile, region: Dict[str, str]) -> Dict[str, str]:
    """Process a single BED region and return the count."""
    return {
        region['name']: str(samfile.count(
            reference=region['chrom'],
            start=int(region['start']),
            end=int(region['end'])
        ))
    }

def extract_counts(datadir: Path, full_BED: Path, sample_id: str, db: Dict) -> None:
    """
    Function that reads the BAM file and extracts the read count for each window

    :param datadir: Location of the BAM files
    :param full_BED: BED file to guide the counts
    :param sample_id: ID of the sample
    :param db: Database for logging

    :type datadir: Path
    :type full_BED: Path
    :type sample_id: str
    :type db: Dict

    :return: None
    """
    console.log(f"Using BED file {full_BED}")
    log_to_db(db, f"Using BED file {full_BED}", "INFO", "count2", sample_id, datadir.name)

    bam_file = datadir / "BAM" / sample_id / "BAM" / f"{sample_id}.bam"
    cnv_file = datadir / "BAM" / sample_id / f"{sample_id}.cnv"

    if not cnv_file.exists():
        try:
            console.log(f"Creating sample's CNV file {cnv_file}")
            log_to_db(db, f"Creating sample's CNV file {cnv_file}", "INFO", "count2", sample_id, datadir.name)

            with open(full_BED, 'r') as bed, \
                 open(cnv_file, 'w') as cnv_out, \
                 pysam.AlignmentFile(bam_file, "rb") as samfile:

                cnv_out.write(f"Location\t{sample_id}\n")

                bed_regions = [{'chrom': line.split('\t')[0],
                                'start': line.split('\t')[1],
                                'end': line.split('\t')[2],
                                'name': line.split('\t')[3]} for line in bed]

                with Progress() as progress:
                    task = progress.add_task("[cyan]Processing regions...", total=len(bed_regions))

                    with ThreadPoolExecutor(max_workers=4) as executor:
                        future_to_region = {executor.submit(process_bed_region, samfile, region): region for region in bed_regions}

                        for future in as_completed(future_to_region):
                            region = future_to_region[future]
                            try:
                                result = future.result()
                                cnv_out.write(f"{region['name']}\t{result[region['name']]}\n")
                            except Exception as exc:
                                console.log(f"Error processing region {region['name']}: {exc}")
                                log_to_db(db, f"Error processing region {region['name']}: {exc}", "ERROR", "count2", sample_id, datadir.name)
                            progress.update(task, advance=1)

            console.log(f"BAM file analyzed, CNV file created for {sample_id}")
            log_to_db(db, f"BAM file analyzed, CNV file created for {sample_id}", "INFO", "count2", sample_id, datadir.name)

            # Save data to ClickHouse
            save_to_clickhouse(cnv_file, sample_id, datadir.name)

        except Exception as e:
            console.log(f"Error: {str(e)}")
            log_to_db(db, f"Error: {str(e)}", "ERROR", "count2", sample_id, datadir.name)
    else:
        console.log(f"CNV file already exists for {sample_id}")
        log_to_db(db, f"CNV file already exists for {sample_id}", "INFO", "count2", sample_id, datadir.name)


def save_to_clickhouse(cnv_file: Path, sample_id: str, run_id: str) -> None:
    client = get_clickhouse_client()

    # Create table if not exists
    client.execute('''
        CREATE TABLE IF NOT EXISTS cnv_data (
            sample_id String,
            run_id String,
            gene String,
            location String,
            count UInt32,
            timestamp DateTime
        ) ENGINE = MergeTree()
        ORDER BY (sample_id, run_id, gene)
    ''')

    # Read CNV file and insert data
    with open(cnv_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            location, count = line.strip().split('\t')
            gene = location.split(':')[0]  # Assuming gene is the first part of the location
            client.execute(
                'INSERT INTO cnv_data (sample_id, run_id, gene, location, count, timestamp) VALUES',
                [(sample_id, run_id, gene, location, int(count), datetime.now())]
            )

    console.log(f"Data for {sample_id} saved to ClickHouse")

@click.command()
@click.option('-d', '--datadir', type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path), required=True, help="Path to the data directory")
@click.option('-b', '--bed', type=str, required=True, help="Name of the BED file to use")
@click.option('-s', '--sample', type=str, required=True, help="Sample ID")
@click.option('-l', '--log', type=click.Path(dir_okay=False, path_type=Path), default='count2_log.json', help="Path to the log file (default: count2_log.json)")
def main(datadir: Path, bed: str, sample: str, log: Path):
    """Extract read counts for genomic regions."""
    # Load environment variables
    load_dotenv()

    # Set up paths and database
    bed_files_dir = os.getenv('BED_FILES')
    if not bed_files_dir:
        raise click.ClickException("BED_FILES environment variable is not set")
    full_BED = Path(bed_files_dir) / bed

    # Validate inputs
    if not full_BED.exists():
        raise click.FileError(str(full_BED), hint="BED file not found")

    # Set up logging
    db = get_sample_db(datadir, sample)

    try:
        click.echo(f"Processing sample {sample} with BED file {full_BED}")
        extract_counts(datadir, full_BED, sample, db)
        click.echo(f"Processing completed for sample {sample}")
    except Exception as e:
        click.echo(f"An error occurred while processing sample {sample}: {str(e)}", err=True)
        raise click.Abort()

if __name__ == "__main__":
    main()