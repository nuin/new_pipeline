# bwa_align.py Documentation

## Overview

`bwa_align.py` is a Python module that performs the alignment/mapping for each pair/sample_id using the BWA software. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `os`
- `subprocess`
- `time`
- `pathlib`
- `typing`
- `rich`

## Functions

The script contains the following function:

- `run_bwa()`: Performs the alignment/mapping for each pair/sample_id. It takes the sample ID, list of FASTQ files for the sample, directory of the run, reference genome, path to the BWA software, and path to the Samtools software as arguments. It returns a flag indicating whether the BAM file exists (1) or was created (0).

## Usage

The `run_bwa()` function is called within the larger NGS pipeline to perform the alignment/mapping for each pair/sample_id. It is not intended to be run as a standalone script.

## Output

The `run_bwa()` function outputs a flag indicating whether the BAM file exists (1) or was created (0). It also logs the start of the BWA processing for the file, whether the BAM file exists or was created, and the BWA and Samtools commands used.

## Example

Here is an example of how the `run_bwa()` function might be used within the larger NGS pipeline:

```python
from src.lib.bwa_align import run_bwa

sample_id = "sample1"
fastq_files = ["sample1_R1.fastq", "sample1_R2.fastq"]
datadir = "/path/to/data"
reference = "/path/to/reference/genome"
bwa = "/path/to/bwa"
samtools = "/path/to/samtools"

bam_exists = run_bwa(sample_id, fastq_files, datadir, reference, bwa, samtools)
```

In this example, the `run_bwa()` function is called with the sample ID, list of FASTQ files for the sample, directory of the run, reference genome, path to the BWA software, and path to the Samtools software. The function returns a flag indicating whether the BAM file exists (1) or was created (0).