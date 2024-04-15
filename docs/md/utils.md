# utils.py Documentation

## Overview

`utils.py` is a Python module that provides utility functions for various operations such as moving BAM files and compiling identity files. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `glob`
- `pathlib`
- `shutil`
- `rich`

## Functions

The script contains the following functions:

- `move_bam()`: Moves the BAM file and its index file from one location to another. It takes the directory where the data is located, the sample ID, and the name of the BAM file as arguments. It returns True if the operation is successful, False otherwise.

- `compile_identity()`: Reads samples' identity files and compiles them in a single file. It takes the directory where the run is located as an argument. It returns True if the operation is successful.

## Usage

The functions in `utils.py` are called within the larger NGS pipeline to perform various utility operations. It is not intended to be run as a standalone script.

## Output

The functions in `utils.py` output various results, including flags indicating the completion of certain functions, and files with compiled identity data.

## Example

Here is an example of how the `move_bam()` and `compile_identity()` functions might be used within the larger NGS pipeline:

```python
from src.lib.utils import move_bam, compile_identity

sample_id = "sample1"
datadir = "/path/to/data"
bam_file = "sample_bam"

move_bam(datadir, sample_id, bam_file)
compile_identity(datadir)
```

In this example, the `move_bam()` function is called with the directory where the data is located, the sample ID, and the name of the BAM file to move the BAM file and its index file from one location to another. The `compile_identity()` function is called with the directory where the run is located to read samples' identity files and compile them in a single file.