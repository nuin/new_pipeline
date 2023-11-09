"""
.. module:: count2.py
    :platform: any
    :synopsis: Module the generate a read count under each exon
.. moduleauthor:: Paulo Nuin, October 2017
"""

from pathlib import Path

import pysam
from rich.console import Console

console = Console()


def extract_counts(datadir, full_BED, sample_id):
    """
    Function that reads the BAM file and extract the read count for each window

    :param datadir: Location of the BAM files
    :param full_BED: BED file to guide the counts
    :param sample_id: ID of the sdample

    :type datadir: string
    :type full_BED: string
    :type sample_id: string

    :return: no return yet
    :todo: change BED file location from hardcoded
    """

    bedfile = open(full_BED).read().splitlines()
    console.log(f"Using BED file {full_BED}")
    bam_file = f"{datadir}/BAM/{sample_id}/BAM/{sample_id}.bam"
    if not Path(f"{datadir}/BAM/{sample_id}/{sample_id}.cnv").exists():
        try:
            console.log(f"Creating samples CNV file {datadir}")
            cnv_out = open(f"{datadir}/BAM/{sample_id}/{sample_id}.cnv", "w")
            cnv_out.write("Location\t")
            cnv_out.write(f"{sample_id}\n")
            console.log(f"Analysing BAM file {sample_id} {datadir}")
            samfile = pysam.AlignmentFile(bam_file, "rb")
            for location in bedfile:
                temp = location.split("\t")
                cnv_out.write(
                    f"{temp[3]}\t{str(samfile.count(reference=temp[0], start=int(temp[1]), end=int(temp[2])))}\n"
                )
            cnv_out.close()

            console.log(f"BAM file analysed, CNV file created {sample_id} {datadir}")
        except Exception as e:
            console.log(str(e))
            console.log(f"BAM file not found {sample_id} {datadir}")
    else:
        console.log(f"CNV file already exists {sample_id} {datadir}")


