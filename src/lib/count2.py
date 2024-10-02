"""
.. module:: count2.py
    :platform: any
    :synopsis: Module to generate read counts under each exon
.. moduleauthor:: Paulo Nuin, October 2017, Updated September 2024
"""

import pysam
from pathlib import Path
from rich.console import Console
from rich.progress import Progress
from .log_api import log_to_api
from .db_logger import log_to_db
from io import StringIO

console = Console()

def extract_counts(datadir: Path, full_BED: Path, sample_id: str, db: dict) -> None:
    """
    Function that reads the BAM file and extracts the read count for each window

    :param datadir: Location of the BAM files
    :param full_BED: BED file to guide the counts
    :param sample_id: ID of the sample
    :param db: Database for logging

    :type datadir: Path
    :type full_BED: Path
    :type sample_id: str
    :type db: dict

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

            with pysam.AlignmentFile(str(bam_file), "rb") as samfile, \
                 open(cnv_file, 'w') as cnv_out, \
                 open(full_BED, 'r') as bed:

                cnv_out.write(f"Location\t{sample_id}\n")
                bed_regions = bed.readlines()

                buffer = StringIO()
                with Progress() as progress:
                    task = progress.add_task("[cyan]Processing regions...", total=len(bed_regions))

                    for location in bed_regions:
                        temp = location.strip().split("\t")
                        try:
                            count = samfile.count(temp[0], int(temp[1]), int(temp[2]))
                            buffer.write(f"{temp[3]}\t{count}\n")
                            if buffer.tell() > 1000000:  # Flush every ~1MB
                                cnv_out.write(buffer.getvalue())
                                buffer.seek(0)
                                buffer.truncate(0)
                        except Exception as e:
                            console.log(f"Error processing region {temp[3]}: {str(e)}")
                            log_to_db(db, f"Error processing region {temp[3]}: {str(e)}", "ERROR", "count2", sample_id, datadir.name)
                            buffer.write(f"{temp[3]}\tERROR\n")
                        progress.update(task, advance=1)

                cnv_out.write(buffer.getvalue())  # Write any remaining data

            console.log(f"BAM file analyzed, CNV file created for {sample_id}")
            log_to_db(db, f"BAM file analyzed, CNV file created for {sample_id}", "INFO", "count2", sample_id, datadir.name)
        except Exception as e:
            console.log(f"Error: {str(e)}")
            log_to_db(db, f"Error: {str(e)}", "ERROR", "count2", sample_id, datadir.name)
    else:
        console.log(f"CNV file already exists for {sample_id}")
        log_to_db(db, f"CNV file already exists for {sample_id}", "INFO", "count2", sample_id, datadir.name)

        
if __name__ == "__main__":
    import click
    import os
    from dotenv import load_dotenv
    from .db_logger import get_sample_db

    load_dotenv()

    @click.command()
    @click.option('-d', '--datadir', type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path), required=True, help="Path to the data directory")
    @click.option('-b', '--bed', type=str, required=True, help="Name of the BED file to use")
    @click.option('-s', '--sample', type=str, required=True, help="Sample ID")
    def main(datadir: Path, bed: str, sample: str):
        """Extract read counts for genomic regions."""
        bed_files_dir = os.getenv('BED_FILES')
        if not bed_files_dir:
            raise click.ClickException("BED_FILES environment variable is not set")
        full_BED = Path(bed_files_dir) / bed

        if not full_BED.exists():
            raise click.FileError(str(full_BED), hint="BED file not found")

        db = get_sample_db(datadir, sample)

        try:
            click.echo(f"Processing sample {sample} with BED file {full_BED}")
            extract_counts(datadir, full_BED, sample, db)
            click.echo(f"Processing completed for sample {sample}")
        except Exception as e:
            click.echo(f"An error occurred while processing sample {sample}: {str(e)}", err=True)
            raise click.Abort()

    main()