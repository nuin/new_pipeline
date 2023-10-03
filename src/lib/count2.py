"""
.. module:: count2.py
    :platform: any
    :synopsis: Module the generate a read count under each exon
.. moduleauthor:: Paulo Nuin, October 2017
"""

import os
import pysam
from pathlib import Path
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
    bam_file = f"{datadir}/BAM/{sample_id}/BAM/{sample_id}.bam"
    # if Path(f"{datadir}/BAM/{sample_id}/{sample_id}.cnv").exists():
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
    # else:
    #     console.log(f"CNV file already exists {sample_id} {datadir}")


if __name__ == "__main__":

    # runs = ['171020_NB551084_0023_AHV2N2AFXX_Cplus_2017_NGS_24']
    # for run in runs:
    #     for sample_id in glob.glob('/Volumes/Jupiter/CancerPlusRuns/' + run + '/BAM/*'):
    #         patient = os.path.basename(sample_id)
    #         extract_counts('/Volumes/Jupiter/CancerPlusRuns/' + run,
    #                        '/opt/BED/Inherited_Cancer_panel_BED_Window_2.bed', patient)
    # extract_counts('/Volumes/Venus/pipeline_validation/CNV/171020_NB551084_0023_AHV2N2AFXX_Cplus_2017_NGS_24',
    # '/opt/BED/Inherited_Cancer_panel_BED_Window_2.bed', '83862_RP')
    extract_counts(
        "/Volumes/Venus/TestDevelopment/PARP/TD_02",
        "/opt/BED/BRCA_amplicon.bed",
        "HD-810_H",
    )
