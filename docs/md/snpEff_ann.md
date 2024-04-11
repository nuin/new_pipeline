# snpEff_ann.py Documentation

## Overview

`snpEff_ann.py` is a Python module that uses SnpEff for annotating and predicting the effects of variants on genes. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following function:

- `snpEff_annotation()`: Calls SnpEff to annotate and predict the effects of variants on genes. It takes the sample ID, directory of the run, reference genome, and SnpEff jar file location as arguments. It returns 'success' if the annotation is successful, 'exists' if the annotated file already exists.

## Usage

The `snpEff_annotation()` function is called within the larger NGS pipeline to annotate and predict the effects of variants on genes. It is not intended to be run as a standalone script.

## Output

The `snpEff_annotation()` function outputs an annotated VCF file. It also logs the start of the annotation process, whether the annotated VCF file exists or was created, and any errors encountered during the process.

## Example

Here is an example of how the `snpEff_annotation()` function might be used within the larger NGS pipeline:

```python
from src.lib.snpEff_ann import snpEff_annotation

sample_id = "sample1"
datadir = "/path/to/data"
reference = "/path/to/reference"
snpEff = "/path/to/snpEff.jar"

result = snpEff_annotation(datadir, sample_id, reference, snpEff)
```

In this example, the `snpEff_annotation()` function is called with the sample ID, directory of the run, reference genome, and SnpEff jar file location. The function returns 'success' if the annotation is successful, 'exists' if the annotated file already exists.