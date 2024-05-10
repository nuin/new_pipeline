# Project Documentation

## Overview

This project is a Next-Generation Sequencing (NGS) pipeline implemented in Python. It performs various operations such as FASTQ alignments, variant calling, and quality control metrics generation. The pipeline is designed to be modular and flexible, allowing for the addition of new steps and modifications of existing ones.

## Dependencies

The project depends on several Python libraries and external tools, including:

- Python libraries: `os`, `subprocess`, `pathlib`, `rich`, `click`, `yaml`, `dotenv`
- External tools: BWA, SAMtools, Picard, GATK, GATK3, Freebayes, Octopus, SnpEff

## Modules

The project is organized into several Python modules, each responsible for a specific part of the pipeline:

- `align.py`: Performs FASTQ alignments using the BWA aligner.
- `utils.py`: Provides utility functions for various operations.
- `variants_freebayes.py`: Generates variants by calling Freebayes.
- `variants_GATK.py`: Generates variants by calling the Genome Analysis Toolkit (GATK).
- `variants_GATK3.py`: Generates variants by calling the Genome Analysis Toolkit version 3 (GATK3).
- `variants_octopus.py`: Generates variants by calling the Octopus variant caller.
- Other modules for specific tasks such as duplicate removal, recalibration, annotation, coverage calculation, and more.

## Usage

The main entry point of the pipeline is the `pipeline.py` script. It uses the Click library to provide a command-line interface. Here is an example of how to run the pipeline:

```bash
python pipeline.py -c config.yaml -d /path/to/data -p panel_name -s sample1 sample2
```

In this example, `config.yaml` is the configuration file for the run, `/path/to/data` is the directory with the location of the FASTQ files, `panel_name` is the panel to be used, and `sample1` and `sample2` are the samples to be analysed.

## Output

The pipeline generates various outputs, including aligned BAM files, variant call files (VCF), and quality control metrics. The outputs are organized in a directory structure with separate directories for each sample and each type of output.

## Configuration

The pipeline is configured through a YAML configuration file. This file specifies the reference genome, the panel to be used, and the paths to the external tools.

## Contributing

Contributions to the project are welcome. Please follow the standard GitHub workflow for contributing to a project.

## License

The project is licensed under the MIT license.