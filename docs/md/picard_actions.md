# picard_actions.py Documentation

## Overview

`picard_actions.py` is a Python module that calls Picard to sort Freebayes VCF. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following function:

- `picard_sort()`: Calls Picard for Freebayes VCF sorting. It takes the sample ID, directory of the run, reference file used in the original alignment, and Picard jar file location as arguments. It returns 'success' if the VCF files are successfully sorted, 'exists' if the sorted file already exists.

## Usage

The `picard_sort()` function is called within the larger NGS pipeline to sort Freebayes VCF. It is not intended to be run as a standalone script.

## Output

The `picard_sort()` function outputs a sorted VCF file. It also logs the start of the VCF sorting process, whether the sorted VCF file exists or was created, and any errors encountered during the process.

## Example

Here is an example of how the `picard_sort()` function might be used within the larger NGS pipeline:

```python
from src.lib.picard_actions import picard_sort

sample_id = "sample1"
datadir = "/path/to/data"
reference = "/path/to/reference"
picard = "/path/to/picard.jar"

result = picard_sort(datadir, sample_id, reference, picard)
```

In this example, the `picard_sort()` function is called with the sample ID, directory of the run, reference file used in the original alignment, and Picard jar file location. The function returns 'success' if the VCF files are successfully sorted, 'exists' if the sorted file already exists.