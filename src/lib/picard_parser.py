"""
.. module:: picard_parser
    :platform: any
    :synopsis: Small module that returns the content of a picard result file, parsed
.. moduleauthor:: Paulo Nuin, June 2016

"""

import sys


def main_parser(picard_output):
    """
    Simple function that parses usual picard output files

    :param picard_output: location of the picard file

    :type picard_output: string

    :return: file contents

    :todo: return error
    """

    picard_file = open(picard_output).readlines()

    datalines = filter(lambda x: not x.startswith("#"), picard_file)
    datalines = filter(lambda x: x != "\n", datalines)
    datalines2 = [x.strip() for x in datalines]

    return datalines2


if __name__ == "__main__":

    main_parser(sys.argv[1])[1]
