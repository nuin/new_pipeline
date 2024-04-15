# count2.py Documentation

## Overview

`count2.py` is a Python module that generates a read count under each exon in Next-Generation Sequencing (NGS) data. It is part of a larger NGS pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `pathlib`
- `pysam`
- `rich`

## Functions

The script contains the following function:

- `extract_counts()`: Reads the BAM file and extracts the read count for each window. It takes the directory of the run, BED file to guide the counts, and the sample ID as arguments. It does not return anything.

## Usage

The `extract_counts()` function is called within the larger NGS pipeline to generate a read count under each exon. It is not intended to be run as a standalone script.

## Output

The `extract_counts()` function outputs a CNV file with the read count for each window. It also logs the start of the processing for the file, whether the CNV file exists or was created, and any errors encountered during the process.

## Example

Here is an example of how the `extract_counts()` function might be used within the larger NGS pipeline:

```python
from src.lib.count2 import extract_counts

datadir = "/path/to/data"
full_BED = "/path/to/bedfile.bed"
sample_id = "sample1"

extract_counts(datadir, full_BED, sample_id)
```

In this example, the `extract_counts()` function is called with the directory of the run, BED file to guide the counts, and the sample ID. The function does not return anything.