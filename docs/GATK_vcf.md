# GATK_vcf.py Documentation

## Overview

`GATK_vcf.py` is a Python module that calls GATK to compare VCF files and performs post-analysis of these variants. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following function:

- `vcf_comparison()`: Merges the available VCFs in the sample VCF datadir. It takes the sample ID, directory of the run, reference file used in the original alignment, and GATK jar file location as arguments. It returns 'success' if the VCF files are successfully merged, 'exists' if the merged file already exists.

## Usage

The `vcf_comparison()` function is called within the larger NGS pipeline to compare VCF files and performs post-analysis of these variants. It is not intended to be run as a standalone script.

## Output

The `vcf_comparison()` function outputs a merged VCF file. It also logs the start of the VCF merging process, whether the merged VCF file exists or was created, and any errors encountered during the process.

## Example

Here is an example of how the `vcf_comparison()` function might be used within the larger NGS pipeline:

```python
from src.lib.GATK_vcf import vcf_comparison

sample_id = "sample1"
datadir = "/path/to/data"
reference = "/path/to/reference"
gatk = "/path/to/gatk.jar"

result = vcf_comparison(datadir, sample_id, reference, gatk)
```

In this example, the `vcf_comparison()` function is called with the sample ID, directory of the run, reference file used in the original alignment, and GATK jar file location. The function returns 'success' if the VCF files are successfully merged, 'exists' if the merged file already exists.