"""
.. module:: variants_table
    :platform: any
    :synopsis: Module that saves a parsed variant table from a series of VCFs
.. moduleauthor:: Paulo Nuin, December 2017
"""

from collections import defaultdict
import vcf
import logging
import yaml
import os
import graypy
from .vcf_parser import parse_annotation

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
graypy_handler = graypy.GELFUDPHandler('192.168.1.169', 12201)
logger.addHandler(graypy_handler)


def get_code(sample_id):
    """

    :param sample_id:
    :return:
    """

    return sample_id[-2:]


def extract_info(sample_id, directory, transcript_location, run_id):
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

    logger.info('Processing table for ' + sample_id)

    transcript_file = open(transcript_location).read().splitlines()

    all_samples = yaml.load(open(directory + '/' + os.path.basename(directory) + '.yaml').read(), Loader=yaml.FullLoader)['Genes']

    sample_genes = all_samples[sample_id].rstrip('_').split('_')

    logger.info('Sample genes are ' + ', '.join(sample_genes) + ' ' + sample_id)

    transcripts, genes = [], []
    transcript_file = open(transcript_location).read().splitlines()

    variants = defaultdict(dict)

    code = get_code(sample_id)
    argument = directory + '/BAM/' + sample_id + '/' + sample_id
    argument_vcf = directory + '/BAM/' + sample_id + '/VCF_' + code + '/' + sample_id

    table_file = open(argument + '_table.txt', 'w')

    annotated_vcf = vcf.Reader(open(argument_vcf + '.hg19_multianno.vcf', 'r'))
    annotated_vcf2 = vcf.Reader(open(argument_vcf + '_merged.ann.vcf', 'r'))

    # get sample
    sample = annotated_vcf.samples[0]

    for record in annotated_vcf2:
        temp_dict = {}
        all_annotation = parse_annotation(record.INFO['ANN'])
        temp_dict['chrom'] = record.CHROM
        temp_dict['pos'] = str(record.POS)
        temp_dict['set'] = record.INFO['set']
        variants[record.CHROM + '_' + str(record.POS)].update(all_annotation)
        variants[record.CHROM + '_' + str(record.POS)].update(temp_dict)

    try:
        for record in annotated_vcf:
            temp_dict = {}
            for info in record.FORMAT.split(':'):
                temp_dict[info] = record.genotype(sample)[info]

            variants[record.CHROM + '_' + str(record.POS)].update(temp_dict)
    except Exception as e:
        logger.error(str(e))

    sorted_variants = sorted(variants.keys())
    table_file.write('CHROM\tPOS\tGene\tTranscript\tc.\tp.\tAllelic ratio\tRead depth\tMAF\tImpact\tCallers\n')
    for variant in sorted_variants:
        try:
            allelic_ratio = str('{0:.2f}'.format(float(variants[variant]['AO']) / float(variants[variant]['DP'])))
            read_depth = str(variants[variant]['DP'])
        except Exception:
            allelic_ratio = 'N/A'
            read_depth = 'N/A'
        if sample_genes == ['ALL']:
            # logger.info('Printing ALL genes ' + variants[variant]['gene'])
            table_file.write(variants[variant]['chrom'] + '\t' + variants[variant]['pos'] + '\t' + variants[variant]['gene'] + '\t')
            table_file.write(variants[variant]['feature_id'] + '\t' + variants[variant]['c_notation'] + '\t' + variants[variant]['p_notation'] + '\t')
            table_file.write(allelic_ratio + '\t' + read_depth + '\t' + 'MAF' + '\t' + variants[variant]['impact'] + '\t' + variants[variant]['set'] + '\n')
        else:
            if variants[variant]['gene'] in sample_genes:
                logger.info('Printing partial ' + variants[variant]['gene'] +  ' ' + sample_id)
                table_file.write(variants[variant]['chrom'] + '\t' + variants[variant]['pos'] + '\t' + variants[variant]['gene'] + '\t')
                table_file.write(variants[variant]['feature_id'] + '\t' + variants[variant]['c_notation'] + '\t' + variants[variant]['p_notation'] + '\t')
                table_file.write(allelic_ratio + '\t' + read_depth + '\t' + 'MAF' + '\t' + variants[variant]['impact'] + '\t' + variants[variant]['set'] + '\n')

    table_file.close()
    return 'success'


if __name__ == '__main__':

    data_directory = '/Volumes/Jupiter/CancerPlusRuns/210731_NB551084_0137_AHMYNKAFY4_Cplus_2021_NGS_15T'
    for sample_id in [
        '21GN-138G00059b_MM_DF',
        '21GN-140G00038b_PS_OS',
        '21GN-143G00001b_DC_HK',
        '21GN-145G00035b_GP_DF',
        '21GN-145G00051b_BF_OS',
        '21GN-145G00053b_BA_OS',
        '21GN-146G00001b_JT_OS',
        '21GN-146G00009b_TC_JA',
        '21GN-146G00020b_MT_1F',
        '21GN-147G00023b_BL_OS',
        '21GN-147G00025b_BI_DF',
        '21GN-147G00036b_GC_DF',
        '21GN-148G00011b_GG_OS',
        '21GN-148G00022b_CT_DF',
        '21GN-148G00033b_KE_HK',
        '21GN-151G00034b_KE_OS',
        '21GN-152G00005b_MK_1F',
        '21GN-160G00036b_KG_1F',
        '21GN-161G00011b_DN_DF',
        '21GN-186G00024b_FM_Z9',
        '21GN-187G00009b_MM_Z9',
        '21GN-187G00029b_KS_DF',
        '21GN-189G00018b_CB_Z9',
        '21GN-189G00027b_CE_DF',
    ]:
        extract_info(sample_id, data_directory, '/opt/bundle/transcripts.txt', 'Test_dataset')
