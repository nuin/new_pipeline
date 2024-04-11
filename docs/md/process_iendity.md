# process_identity.py Documentation

## Overview

`process_identity.py` is a Python module that processes identity files and returns a full table with barcodes and other information. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `glob`
- `os`
- `pathlib`
- `numpy`
- `pandas`
- `rich`

## Functions

The script contains the following functions:

- `read_identity()`: Reads the full identity for each patient and returns a DataFrame with information. It takes the sample ID and run location as arguments.

- `generate_barcode()`: Reads the codes and generates a numeric barcode for each sample. It takes the sample ID, run location, and a DataFrame with sample identity information as arguments.

- `process_identity()`: Processes the identity for each patient and returns a DataFrame with information. It takes the sample ID, run location, and a DataFrame with sample identity information as arguments.

- `barcoding()`: Processes the identity for each patient, generates a numeric barcode for each sample, and saves the results. It takes the sample ID and run location as arguments.

- `find_duplicates()`: Checks for duplicate barcodes in the data directory. It takes the run location as an argument.

- `compile_barcodes()`: Compiles all barcodes from the data directory into a single file. It takes the run location as an argument.

## Usage

The functions in `process_identity.py` are called within the larger NGS pipeline to process identity files and generate barcodes. It is not intended to be run as a standalone script.

## Output

The functions in `process_identity.py` output various results, including DataFrames with sample data, numeric barcodes for each sample, and files with barcodes.

## Example

Here is an example of how the `read_identity()`, `process_identity()`, and `generate_barcode()` functions might be used within the larger NGS pipeline:

```python
from src.lib.process_identity import read_identity, process_identity, generate_barcode

sample_id = "sample1"
datadir = "/path/to/data"

sample_identity = read_identity(sample_id, datadir)
processed_identity = process_identity(sample_id, datadir, sample_identity)
generate_barcode(sample_id, datadir, processed_identity)
```

In this example, the `read_identity()` function is called with the sample ID and run location to read the full identity for the patient. The `process_identity()` function is called with the sample ID, run location, and the DataFrame with sample identity information to process the identity. The `generate_barcode()` function is called with the sample ID, run location, and the DataFrame with processed sample identity information to generate a numeric barcode for the sample.