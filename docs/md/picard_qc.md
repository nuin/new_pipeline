# picard_qc.py Documentation

## Overview

`picard_qc.py` is a Python module that generates nucleotide-based coverage using Picard. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `os`
- `string`
- `subprocess`
- `pathlib`
- `rich`
- `suds`

## Functions

The script contains the following functions:

- `get_coverage()`: Calls Picard to generate nucleotide coverage. It takes the sample ID, directory of the run, reference genome, Picard specific BED file, Picard jar file location, and panel type as arguments. It returns 'success' if the Picard coverage file is successfully created, 'exists' if the file already exists.

- `get_coverage_parp()`: Calls Picard to generate nucleotide coverage for PARP. It takes the sample ID, directory of the run, reference genome, Picard specific BED file, and Picard jar file location as arguments. It returns 'success' if the Picard coverage file is successfully created, 'exists' if the file already exists.

- `get_transcripts()`: Gets the currently used transcripts for HGVS numberConversion. It takes the location of the transcript file as an argument. It returns a dictionary of transcripts.

- `chr_frame()`: Checks chromosome regions for coverage under 25X. It takes the current segment being analysed as an argument. It returns the gene, segment and chromosome if there's a region under 25x, otherwise returns 'empty', 'empty', 'empty'.

- `convert_g_to_c()`: Converts g. notation to c. It takes the chromosome number, position to be converted, and list of transcripts as arguments. It returns the HVGS notation of the position.

- `create_table()`: Creates a text table with the locations under 20x of coverage. It takes a list of regions under 100x, chromosome number, gene symbol, transcript ID, and location of the coverage file as arguments. It does not return anything.

## Usage

The functions in `picard_qc.py` are called within the larger NGS pipeline to generate nucleotide-based coverage using Picard. It is not intended to be run as a standalone script.

## Output

The functions in `picard_qc.py` output various results, including flags indicating the completion of certain functions, and files with Picard coverage.

## Example

Here is an example of how the `get_coverage()` and `get_coverage_parp()` functions might be used within the larger NGS pipeline:

```python
from src.lib.picard_qc import get_coverage, get_coverage_parp

sample_id = "sample1"
datadir = "/path/to/data"
reference = "/path/to/reference"
picard = "/path/to/picard.jar"
bait_file = "/path/to/bait_file"

result_coverage = get_coverage(sample_id, datadir, reference, bait_file, picard)
result_coverage_parp = get_coverage_parp(sample_id, datadir, reference, bait_file, picard)
```

In this example, the `get_coverage()` function is called with the sample ID, directory of the run, reference genome, Picard specific BED file, Picard jar file location, and panel type. The function returns 'success' if the Picard coverage file is successfully created, 'exists' if the file already exists. The `get_coverage_parp()` function is called with the sample ID, directory of the run, reference genome, Picard specific BED file, and Picard jar file location. The function returns 'success' if the Picard coverage file is successfully created, 'exists' if the file already exists.