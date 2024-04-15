# dup_indels.py Documentation

## Overview

`dup_indels.py` is a Python module that performs two calls for Picard to remove duplicate reads and to add header information to the BAM files so they can be analysed in the remainder of the pipeline. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following function:

- `remove_duplicates()`: Runs the duplicate removal step in the analysis using Picard. It takes the sample ID, directory of the run, and Picard jar file location as arguments. It returns 'success' if the process completes successfully, 'exists' if the dedup.bam file already exists, or 'error - process' if the duplicate removal fails.

## Usage

The `remove_duplicates()` function is called within the larger NGS pipeline to remove duplicate reads and to add header information to the BAM files. It is not intended to be run as a standalone script.

## Output

The `remove_duplicates()` function outputs a dedup.bam file with the duplicate reads removed. It also logs the start of the duplicate removal process, whether the dedup.bam file exists or was created, and any errors encountered during the process.

## Example

Here is an example of how the `remove_duplicates()` function might be used within the larger NGS pipeline:

```python
from src.lib.dup_indels import remove_duplicates

sample_id = "sample1"
datadir = "/path/to/data"
picard = "/path/to/picard.jar"

result = remove_duplicates(sample_id, datadir, picard)
```

In this example, the `remove_duplicates()` function is called with the sample ID, directory of the run, and Picard jar file location. The function returns 'success' if the process completes successfully, 'exists' if the dedup.bam file already exists, or 'error - process' if the duplicate removal fails.