# uniformity.py Documentation

## Overview

`uniformity.py` is a Python module that calculates the uniformity of coverage for Next-Generation Sequencing (NGS) data. It is part of a larger NGS pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `numpy`
- `pandas`
- `pathlib`
- `rich`

## Functions

The script contains the following functions:

- `calculate_uniformity()`: Calculates the uniformity of coverage for a given sample. It takes the sample ID, directory of the run, and the coverage file as arguments. It returns a DataFrame with the uniformity metrics.

- `write_uniformity()`: Writes the uniformity metrics to a file. It takes the sample ID, directory of the run, and the DataFrame with the uniformity metrics as arguments. It does not return anything.

## Usage

The functions in `uniformity.py` are called within the larger NGS pipeline to calculate and write the uniformity of coverage. It is not intended to be run as a standalone script.

## Output

The functions in `uniformity.py` output a DataFrame with the uniformity metrics and a file with these metrics.

## Example

Here is an example of how the `calculate_uniformity()` and `write_uniformity()` functions might be used within the larger NGS pipeline:

```python
from src.lib.uniformity import calculate_uniformity, write_uniformity

sample_id = "sample1"
datadir = "/path/to/data"
coverage_file = "/path/to/coverage_file"

uniformity_metrics = calculate_uniformity(sample_id, datadir, coverage_file)
write_uniformity(sample_id, datadir, uniformity_metrics)
```

In this example, the `calculate_uniformity()` function is called with the sample ID, directory of the run, and the coverage file to calculate the uniformity of coverage. The `write_uniformity()` function is called with the sample ID, directory of the run, and the DataFrame with the uniformity metrics to write these metrics to a file.