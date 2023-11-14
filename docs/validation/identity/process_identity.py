"""
.. module:: process_identity
    :platform: any
    :synopsis: TBA
.. moduleauthor:: Paulo Nuin, March 2018

"""

# pylint


import sys
import os
import subprocess
import pandas as pd
import numpy as np
import socket
import glob
# from .kafka_logger import log_to_kafka

my_hostname = socket.gethostbyname(socket.gethostname())


genes = {   29416481:  'ALK_E29',
            215645464: 'BARD1_E4',
            14187449:  'XPC_E16',
            112164561: 'APC_E16',
            55214348:  'EGFR_E4',
            116435768: 'MET_E20',
            116436022: 'MET_E21',
            30999280:  'WRN_E26',
            90967711:  'NBN_E10',
            145742879: 'RECQL4_E3',
            98209594:  'PTCH1_E23',
            133208979: 'POLE_E45',
            32936646:  'BRCA2_E17',
            89858602:  'FANCI_E37',
            89849480:  'FANCA_E16',
            59763347:  'BRIP1_E19'
            }

hgvs = {    29416481:  'c.4472A>G',
            215645464: 'c.1134G>C',
            14187449:  'c.2815C>A',
            112164561: 'c.1635G>A',
            55214348:  'c.474C>T',
            116435768: 'c.3912C>T',
            116436022: 'c.4071G>A',
            30999280:  'c.3222G>T',
            90967711:  'c.1197T>C',
            145742879: 'c.132A>G',
            98209594:  'c.3944C>T',
            133208979: 'c.6252A>G',
            32936646:  'c.7806-14T>C',
            89858602:  'c.3906T>C',
            89849480:  'c.1501G>A',
            59763347:  'c.2755T>C'
            }


codes = { 'ALK_E29000T': 1, 'ALK_E290C0T': 2, 'ALK_E290C00': 3, 'ALK_E29': 5,
          'BARD1_E40C00':1, 'BARD1_E40CG0':2, 'BARD1_E400G0':3, 'BARD1_E4': 5,
          'XPC_E1600G0': 1, 'XPC_E1600GT': 2, 'XPC_E16000T': 3, 'XPC_E16': 5,
          'APC_E1600G0': 1, 'APC_E16A0G0': 2, 'APC_E16A000': 3, 'APC_E16': 5,
          'EGFR_E40C00': 1, 'EGFR_E40C0T': 2, 'EGFR_E4000T': 3, 'EGFR_E4': 5,
          'MET_E200C00': 1, 'MET_E200C0T': 2, 'MET_E20000T': 3, 'MET_E20': 5,
          'MET_E2100G0': 1, 'MET_E21A0G0': 2, 'MET_E21A000': 3, 'MET_E21': 5,
          'WRN_E2600G0': 1, 'WRN_E2600GT': 2, 'WRN_E26000T': 3, 'WRN_E26': 5,
          'NBN_E10A000': 1, 'NBN_E10A0G0': 2, 'NBN_E1000G0': 3, 'NBN_E10': 5,
          'RECQL4_E3000T': 1, 'RECQL4_E30C0T': 2, 'RECQL4_E30C00': 3, 'RECQL4_E3': 5,
          'PTCH1_E2300G0': 1, 'PTCH1_E23A0G0': 2, 'PTCH1_E23A000': 3, 'PTCH1_E23': 5,
          'POLE_E45000T':  1, 'POLE_E450C0T':  2, 'POLE_E450C00':  3, 'POLE_E45':  5,
          'BRCA2_E17000T': 1,'BRCA2_E170C0T': 2,'BRCA2_E170C00':  3,'BRCA2_E17':  5,
          'FANCI_E37000T': 1, 'FANCI_E370C0T': 2, 'FANCI_E370C00': 3, 'FANCI_E37': 5,
          'FANCA_E160C00': 1,'FANCA_E160C0T':  2,'FANCA_E16000T':  3,'FANCA_E16':  5,
          'BRIP1_E19A000': 1,'BRIP1_E19A0G0':  2,'BRIP1_E1900G0':  3,'BRIP1_E19':  5
          }

def read_identity(patient_id, directory):

    patient_identity =  pd.read_csv(directory + '/BAM/' + patient_id + '/identity.txt', sep='\t', names = ['Chromosome', 'Position', 'A', 'C', 'G', 'T'])


    return patient_identity


def generate_barcode(patient_id, directory, patient_identity):

    barcode = ''
    barcode_file = open(directory + '/BAM/' + patient_id + '/' + patient_id + '.barcode.txt', 'w')
    for i in [0, 3, 1, 12, 15, 4, 14, 13, 5, 6, 8, 11, 10, 9, 7, 2]:
        barcode += str(codes[patient_identity['Full code'][i]])

    barcode_file.write(barcode)
    print('INFO\t\tBarcode file generated')

def process_identity(patient_id, directory, patient_identity):

    patient_identity['ID'] = patient_identity['Position'].map(genes)
    patient_identity['HGVS'] = patient_identity['Position'].map(hgvs)
    patient_identity['Total reads'] = patient_identity[['A', 'C', 'G', 'T']].sum(axis=1)
    patient_identity['pcA'] = patient_identity['A']/patient_identity['Total reads']
    patient_identity['pcC'] = patient_identity['C']/patient_identity['Total reads']
    patient_identity['pcG'] = patient_identity['G']/patient_identity['Total reads']
    patient_identity['pcT'] = patient_identity['T']/patient_identity['Total reads']
    patient_identity['codeA'] = np.where(patient_identity['pcA'] > 0.1, 'A', '0')
    patient_identity['codeC'] = np.where(patient_identity['pcC'] > 0.1, 'C', '0')
    patient_identity['codeG'] = np.where(patient_identity['pcG'] > 0.1, 'G', '0')
    patient_identity['codeT'] = np.where(patient_identity['pcT'] > 0.1, 'T', '0')
    patient_identity['Full code'] = patient_identity['ID'] + patient_identity['codeA'] + patient_identity['codeC'] + patient_identity['codeG'] + patient_identity['codeT']

    patient_identity.to_csv(directory + '/BAM/' + patient_id + '/' + patient_id + '.identity_full.txt', sep='\t', index=False)
    print('INFO\t\tFull identity file generated')

    return patient_identity

def barcoding(patient_id, directory):

    patient_identity = read_identity(patient_id, directory)
    process_identity(patient_id, directory, patient_identity)
    generate_barcode(patient_id, directory, patient_identity)

def find_duplicates(directory):

  barcodes = open(directory + '/barcodes.txt').read().splitlines()

  unique = set([])
  samples = []
  for line in barcodes:
    temp = line.split('\t')
    samples.append(temp[0])
    unique.add(temp[-1])

  duplicate_dict = {}

  if len(unique) == len(barcodes):
    for sample in samples:
        duplicate_dict[sample] = False
    return duplicate_dict
  else:
    return 'needs implementation'


def compile_barcodes(directory):

  barcodes = open(directory + '/barcodes.txt', 'w')
  barcode_dict = {}
  for sample in glob.glob(directory + '/BAM/*'):
    patient_id = os.path.basename(sample)
    barcode_file = open(directory + '/BAM/' + patient_id + '/' + patient_id + '.barcode.txt').read().splitlines()
    barcode_dict[patient_id] = barcode_file[0]
    barcodes.write(patient_id + '\t' + barcode_file[0] + '\n')
  print('INFO\t\tBarcodes file created')

  barcodes.close()

  # print(find_duplicates(directory))
  return 'success'



if __name__ == '__main__':

    patient_id = '17-341-020232B_WJ_Z9'
    # directory = '/Volumes/Jupiter/CancerPlusRuns/180223_NB551084_0031_AHW3YVAFXX_Cplus_2018_NGS_06/'
    directory = '/Volumes/Jupiter/CancerPlusRuns/180309_NB551084_0032_AHW3WYAFXX_Cplus_2018_NGS_08/'

    # barcoding(patient_id, directory)
    compile_barcodes(directory)
