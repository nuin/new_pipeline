# cnv.py Documentation

## Overview

`cnv.py` is a Python module that calculates and generates normalized read information for Copy Number Variations (CNVs) in Next-Generation Sequencing (NGS) data. It is part of a larger NGS pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `glob`
- `os`
- `sys`
- `pathlib`
- `typing`
- `pandas`
- `requests`
- `yaml`
- `rich`

## Functions

The script contains several functions, including:

- `print_full()`: Prints the full DataFrame without truncation.
- `check_file_size()`: Checks if the number of lines in the CNV file matches the expected number of lines.
- `compile_samples()`: Checks if all .cnv files were generated and extracts the information in each file, returning a pandas DataFrame.
- `cnv_calculation()`: Generates the CNV intra and inter-sample normalization.
- `split_genders()`: Returns a list of males and females for X-linked and other gender related checks.
- `get_xlinked()`: Sends a GET request to a specific URL and returns the JSON response.
- `cnvs_calculation_xlinked()`: Calculates CNVs for x-linked genes.

## Usage

The functions in `cnv.py` are called within the larger NGS pipeline to calculate and generate normalized read information for CNVs. It is not intended to be run as a standalone script.

## Output

The functions in `cnv.py` output various results, including flags indicating the completion of certain functions, pandas DataFrames with CNV data, and files with CNV data and calculations.

## Example

Here is an example of how the `cnv_calculation()` function might be used within the larger NGS pipeline:

```python
from src.lib.cnv import cnv_calculation

datadir = "/path/to/data"
cnvs = pd.DataFrame({"sample1": [1, 2, 3], "sample2": [4, 5, 6]})
yaml_file = "/path/to/config.yaml"

result = cnv_calculation(datadir, cnvs, yaml_file)
```

In this example, the `cnv_calculation()` function is called with the directory of the run, a DataFrame with raw coverage data for each window, and the run configuration file. The function returns a string indicating the completion of the function.