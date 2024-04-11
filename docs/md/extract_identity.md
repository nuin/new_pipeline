# extract_identity.py Documentation

## Overview

`extract_identity.py` is a Python module that extracts information on the nucleotide in the 16 locations that determine identity. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `os`
- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following functions:

- `mpileup()`: Calls samtools to generate the pileup file that can be searched for the nucleotide state in the locations of interest. It takes the sample ID, directory of the run, identity file name with defined location, and Samtools executable location as arguments. It returns 'success' if the pileup file is successfully generated, 'exists' if the file already exists.

- `get_nucleotides()`: Returns nucleotide counts for the region. It takes the reads under a defined location as an argument.

- `create_identity_table()`: Creates the identity table in the sample data directory. It takes the sample ID and location where the table will be saved as arguments.

## Usage

The functions in `extract_identity.py` are called within the larger NGS pipeline to extract information on the nucleotide in the 16 locations that determine identity. It is not intended to be run as a standalone script.

## Output

The functions in `extract_identity.py` output various results, including flags indicating the completion of certain functions, and files with identity information.

## Example

Here is an example of how the `mpileup()` and `create_identity_table()` functions might be used within the larger NGS pipeline:

```python
from src.lib.extract_identity import mpileup, create_identity_table

sample_id = "sample1"
datadir = "/path/to/data"
identity = "/path/to/identity"
samtools = "/path/to/samtools"

result = mpileup(sample_id, datadir, identity, samtools)
create_identity_table(sample_id, datadir)
```

In this example, the `mpileup()` function is called with the sample ID, directory of the run, identity file name with defined location, and Samtools executable location. The function returns 'success' if the pileup file is successfully generated, 'exists' if the file already exists. The `create_identity_table()` function is called with the sample ID and location where the table will be saved.