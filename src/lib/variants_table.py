"""
.. module:: variants_table
    :platform: any
    :synopsis: Module that saves a parsed variant table from a series of VCFs
.. moduleauthor:: Paulo Nuin, December 2017
"""

from collections import defaultdict
import vcf
import yaml
import os

# from .vcf_parser import parse_annotation

from pathlib import Path

from rich.console import Console

console = Console()


def get_code(sample_id):
    """

    :param sample_id:
    :return:
    """

    return sample_id[-2:]


def extract_info(sample_id, datadir, configuration):
    """
    Function that extract data from a annotated (snpEff) VCF file and parse into a table for eacy reading
    and comnparison

    :param sample_id: ID of the sample
    :param directory: run location
    :param transcript_location: location of the transcript file (hardcoded)
    :param run_id: ID of the run

    :type sample_id: string
    :type directory: string
    :type transcript_location: string
    :type run_id: string

    """

    console.log(f"Processing table for {sample_id}")

    # print(configuration)

    genes = configuration["Genes"][sample_id].split("_")

    print(genes)

    console.log(f"Sample genes are {', '.join(genes)} {sample_id}")
    #
    variants = defaultdict(dict)
    #
    argument = directory + "/BAM/" + sample_id + "/" + sample_id
    # argument_vcf = directory + "/BAM/" + sample_id + "/VCF_" + code + "/" + sample_id

    vcf_dir = f"{datadir}/BAM/{sample_id}/VCF}"
    table_file = open(f"{dadatadir}/BAM/{sample_id}/{sample_id}_table.txt", "w")
    annotated_vcf = vcf.Reader(open(f"{vcf_dir}/{sample_id}_merged.ann.vcf", "r"))


    # # get sample
    # sample = annotated_vcf.samples[0]
    #
    # for record in annotated_vcf2:
    #     temp_dict = {}
    #     all_annotation = parse_annotation(record.INFO["ANN"])
    #     temp_dict["chrom"] = record.CHROM
    #     temp_dict["pos"] = str(record.POS)
    #     temp_dict["set"] = record.INFO["set"]
    #     variants[record.CHROM + "_" + str(record.POS)].update(all_annotation)
    #     variants[record.CHROM + "_" + str(record.POS)].update(temp_dict)
    #
    # try:
    #     for record in annotated_vcf:
    #         temp_dict = {}
    #         for info in record.FORMAT.split(":"):
    #             temp_dict[info] = record.genotype(sample)[info]
    #
    #         variants[record.CHROM + "_" + str(record.POS)].update(temp_dict)
    # except Exception as e:
    #     logger.error(str(e))
    #
    # sorted_variants = sorted(variants.keys())
    # table_file.write(
    #     "CHROM\tPOS\tGene\tTranscript\tc.\tp.\tAllelic ratio\tRead depth\tMAF\tImpact\tCallers\n"
    # )
    # for variant in sorted_variants:
    #     try:
    #         allelic_ratio = str(
    #             "{0:.2f}".format(
    #                 float(variants[variant]["AO"]) / float(variants[variant]["DP"])
    #             )
    #         )
    #         read_depth = str(variants[variant]["DP"])
    #     except Exception:
    #         allelic_ratio = "N/A"
    #         read_depth = "N/A"
    #     if sample_genes == ["ALL"]:
    #         # logger.info('Printing ALL genes ' + variants[variant]['gene'])
    #         table_file.write(
    #             variants[variant]["chrom"]
    #             + "\t"
    #             + variants[variant]["pos"]
    #             + "\t"
    #             + variants[variant]["gene"]
    #             + "\t"
    #         )
    #         table_file.write(
    #             variants[variant]["feature_id"]
    #             + "\t"
    #             + variants[variant]["c_notation"]
    #             + "\t"
    #             + variants[variant]["p_notation"]
    #             + "\t"
    #         )
    #         table_file.write(
    #             allelic_ratio
    #             + "\t"
    #             + read_depth
    #             + "\t"
    #             + "MAF"
    #             + "\t"
    #             + variants[variant]["impact"]
    #             + "\t"
    #             + variants[variant]["set"]
    #             + "\n"
    #         )
    #     else:
    #         if variants[variant]["gene"] in sample_genes:
    #             logger.info(
    #                 "Printing partial " + variants[variant]["gene"] + " " + sample_id
    #             )
    #             table_file.write(
    #                 variants[variant]["chrom"]
    #                 + "\t"
    #                 + variants[variant]["pos"]
    #                 + "\t"
    #                 + variants[variant]["gene"]
    #                 + "\t"
    #             )
    #             table_file.write(
    #                 variants[variant]["feature_id"]
    #                 + "\t"
    #                 + variants[variant]["c_notation"]
    #                 + "\t"
    #                 + variants[variant]["p_notation"]
    #                 + "\t"
    #             )
    #             table_file.write(
    #                 allelic_ratio
    #                 + "\t"
    #                 + read_depth
    #                 + "\t"
    #                 + "MAF"
    #                 + "\t"
    #                 + variants[variant]["impact"]
    #                 + "\t"
    #                 + variants[variant]["set"]
    #                 + "\n"
    #             )
    #
    # table_file.close()
    return "success"


if __name__ == "__main__":

    data_directory = "/Volumes/Jupiter/CardioRuns/230921_NB551084_0269_AHJ3H3AFX5_Cardiac_2023_NGS_41"
    for sample_id in [
        "23GN-230G00029",
    ]:
        extract_info(
            sample_id, data_directory, "/opt/bundle/transcripts.txt", "Test_dataset"
        )
