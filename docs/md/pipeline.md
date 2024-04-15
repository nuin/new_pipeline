# Pipeline.py Documentation

## Overview

`pipeline.py` is a Python script that orchestrates a Next-Generation Sequencing (NGS) pipeline. It performs FASTQ alignments and calls all other steps in the process. The script is designed to be run from the command line and takes several command-line arguments.

## Dependencies

The script depends on several Python libraries, including:

- `os`
- `collections`
- `pathlib`
- `typing`
- `click`
- `yaml`
- `dotenv`
- `rich`

It also depends on several external tools, the paths to which are specified in a `.env` file:

- BWA
- SAMTOOLS
- PICARD
- GATK
- GATK3
- FREEBAYES
- SNPEFF
- OCTOPUS

## Functions

The script contains several functions, including:

- `split_string()`: Splits a string into a list of items.
- `find_fastq()`: Finds all the FASTQ files in the given directory.
- `get_ids()`: Extracts the sample IDs from the list of FASTQ files.
- `create_directories()`: Creates necessary directories for each sample in the given directory.
- `align_files()`: Aligns the FASTQ files for each sample using the BWA aligner.
- `process_dir()`: Main function that orchestrates the entire pipeline of processing the data.
- `analyse_pairs()`: Orchestrates the entire pipeline of processing the data for each sample.
- `generate_analysis()`: Main function of the script, find input files, start alignment and processing.
- `run_analysis()`: Starts the analysis process by calling the `generate_analysis` function.

## Usage

The script is run from the command line with the following syntax:

```bash
python pipeline.py -c <configuration_file> -d <datadir> -p <panel> -s <sample1> -s <sample2> ...
```

Where:

- `<configuration_file>` is the YAML configuration file for the run.
- `<datadir>` is the directory of the run.
- `<panel>` is the panel to be used.
- `<sample1>`, `<sample2>`, etc. are the samples to be analysed. If no samples are provided, all samples in the run will be analysed.

## Output

The script outputs a dictionary with the success/failure of each step for each sample.