"""
.. module:: extract_identity
    :platform: Any
    :synopsis: This module extracts information on the nucleotide in the 16 locations that determine identity
.. moduleauthor:: Paulo Nuin, August 2016
"""

# :pylint: 7.0

from __future__ import print_function
import os
import subprocess
import socket
from .kafka_logger import log_to_kafka

my_hostname = socket.gethostbyname(socket.gethostname())

# GENES = {'ALK': ('chr2', 29416481),
#          'APC': ('chr5', 112164561),
#          'BARD1': ('chr2', 215645464),
#          'BRCA2': ('chr13', 32936646),
#          'BRIP1': ('chr17', 59763347),
#          'EGFR': ('chr7', 55214348),
#          'FANCA': ('chr16', 89849480),
#          'FANCI': ('chr15', 89858602),
#          'MET': ('chr7', 116435768),
#          'MET2': ('chr7', 116436022),
#          'NBN': ('chr8', 90967711),
#          'POLE': ('chr12', 133208979),
#          'PTCH1': ('chr9', 98209594),
#          'RECQL4': ('chr8', 145742879),
#          'WRN': ('chr8', 30999280),
#          'XPC': ('chr3', 14187449)}


def mpileup(patient_id, directory, identity, samtools):
    """
    Function that call sammtool to generate the pileup
    file that can be searched for the nucleotide state in
    the locations of interest

    :param patient_id: ID of the patient/sample being analysed
    :param directory: Location of the BAM files
    :param identity: Indentity file name with defined location
    :param samtools: Samtools executable location

    :type patient_id: string
    :type directory: string
    :type identity: string
    :type samtools: string

    :return: returns success or exists

    :todo: return error
    """

    argument2 = directory + '/BAM/' + patient_id

    if os.path.isdir(directory + '/BAM/' + patient_id + '/BAM'):
        argument = directory + '/BAM/' + patient_id + '/BAM/' + patient_id
    else:
        argument = directory + '/BAM/' + patient_id + '/' + patient_id

    if os.path.isfile(argument2 + '/identity.mpileup'):
        log_to_kafka('Identity mpileup exists ' + patient_id + ' ' + my_hostname)
        print('INFO\t\tIdentity mpileup file exists')
        return 'exists'
    else:
        log_to_kafka('Starting mpileup process for sample ' + patient_id + ' ' + my_hostname)
        print('INFO\t\tStarting mpileup process for identity file')
        mpileup_string = '%s mpileup -l %s %s.recal_reads.bam > %s/identity.mpileup' % (samtools, identity, argument, argument2)
        proc = subprocess.Popen(mpileup_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
        log_to_kafka('Mpileup file generation completed ' + patient_id + ' ' + my_hostname)
        print('INFO\t\tmpileup generation completed')
        return 'success'


def get_nucleotides(pileup):
    """
    Function that returns nucleotide counts for the region

    :param pileup: reads under a defined location

    :type pileup: string

    :returns: returns the number of A, C, G and T

    """

    reads = ''.join(e for e in pileup if e.isalnum()).upper()

    return reads.count('A'), reads.count('C'), reads.count('G'), reads.count('T')


def create_identity_table(patient_id, directory):
    """
    Function that creates the identity table in the sample
    directory

    :param patient_id: patient id to be analyses
    :param directory: directory where the table will be saved

    :type patient_id: string
    :type directory: string

    """

    mpileup = open(directory + '/BAM/' + patient_id +'/identity.mpileup').read().splitlines()

    if os.path.isfile(directory + '/BAM/' + patient_id + '/identity.txt'):
        print('INFO\t\tIdentity file exists')
    else:
        print('INFO\t\tStarting creation of identity file')
        identity = open(directory + '/BAM/' + patient_id + '/identity.txt', 'w')
        for line in mpileup:
            temp = line.split()
            # print(line)
            nucleotides = get_nucleotides(temp[4])
            identity.write(temp[0] + '\t' + temp[1] +
                           '\t' + str(nucleotides[0]) + '\t')
            identity.write(str(nucleotides[1]) +
                           '\t' + str(nucleotides[2]) + '\t')
            identity.write(str(nucleotides[3]) + '\n')
        print('INFO\t\tIdentity file creation complete')


if __name__ == '__main__':

    directory = '/Users/nuin/New_projects/Data/Illumina/'
    patient_id = 'NA12877_1'
    identity_file = 'identity.txt'
    samtools = 'samtools'
    mpileup(patient_id, directory, identity_file, samtools)
    create_identity_table(patient_id, directory)
