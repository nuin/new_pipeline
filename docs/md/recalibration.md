# recalibration.py Documentation

## Overview

`recalibration.py` is a Python module that performs base quality recalibration as part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following functions:

- `base_recal1()`: Performs the first step of base recalibration, creating a recalibration data table. It takes the sample ID, run location, BED file with regions to be analysed, VCF file of known regions of variants, reference file used in the original alignment, and GATK jar file location as arguments. It returns 'success' if the operation is successful, 'exists' if the recalibration data table already exists.

- `recalibrate()`: Performs the third step of the base recalibration process, generating the final BAM file. It takes the sample ID, run location, reference file used in the original alignment, and GATK jar file location as arguments. It returns 'success' if the operation is successful, 'exists' if the recalibration data table already exists.

## Usage

The functions in `recalibration.py` are called within the larger NGS pipeline to perform base quality recalibration. It is not intended to be run as a standalone script.

## Output

The functions in `recalibration.py` output various results, including flags indicating the completion of certain functions, and files with recalibration data.

## Example

Here is an example of how the `base_recal1()` and `recalibrate()` functions might be used within the larger NGS pipeline:

```python
from src.lib.recalibration import base_recal1, recalibrate

sample_id = "sample1"
datadir = "/path/to/data"
bed_file = "/path/to/bed_file"
vcf_file = "/path/to/vcf_file"
reference = "/path/to/reference"
gatk = "/path/to/gatk.jar"

base_recal1(datadir, sample_id, bed_file, vcf_file, reference, gatk)
recalibrate(datadir, sample_id, reference, gatk)
```

In this example, the `base_recal1()` function is called with the sample ID, run location, BED file with regions to be analysed, VCF file of known regions of variants, reference file used in the original alignment, and GATK jar file location. The function returns 'success' if the operation is successful, 'exists' if the recalibration data table already exists. The `recalibrate()` function is called with the sample ID, run location, reference file used in the original alignment, and GATK jar file location. The function returns 'success' if the operation is successful, 'exists' if the recalibration data table already exists.