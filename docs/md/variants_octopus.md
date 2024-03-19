# variants_octopus.py Documentation

## Overview

`variants_octopus.py` is a Python module that generates variants by calling the Octopus variant caller. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `os`
- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following functions:

- `octopus_caller()`: Calls Octopus to generate a VCF file. It takes the directory where the data is located, the sample ID, the reference file used in the original alignment, the BED file with regions to be analysed, and the location of Octopus executable as arguments. It returns 'success' if the operation is successful, 'exists' if the VCF file already exists, 'error' if there is an error.

- `edit_octopus_vcf()`: Removes extra lines in Octopus generated VCF to allow proper sorting. It takes the sample ID and the directory where the data is located as arguments. It returns 'exists' if the VCF file already exists, 'success' if the operation is successful.

## Usage

The functions in `variants_octopus.py` are called within the larger NGS pipeline to generate variants and edit the generated VCF file. It is not intended to be run as a standalone script.

## Output

The functions in `variants_octopus.py` output various results, including flags indicating the completion of certain functions, and VCF files with variants.

## Example

Here is an example of how the `octopus_caller()` and `edit_octopus_vcf()` functions might be used within the larger NGS pipeline:

```python
from src.lib.variants_octopus import octopus_caller, edit_octopus_vcf

sample_id = "sample1"
datadir = "/path/to/data"
reference = "/path/to/reference"
bed_file = "/path/to/bed_file"
octopus = "/path/to/octopus"

octopus_caller(datadir, sample_id, reference, bed_file, octopus)
edit_octopus_vcf(sample_id, datadir)
```

In this example, the `octopus_caller()` function is called with the directory where the data is located, the sample ID, the reference file used in the original alignment, the BED file with regions to be analysed, and the location of Octopus executable to generate a VCF file. The `edit_octopus_vcf()` function is called with the sample ID and the directory where the data is located to remove extra lines in Octopus generated VCF to allow proper sorting.