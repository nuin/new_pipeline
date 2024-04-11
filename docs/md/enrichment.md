# enrichment.py Documentation

## Overview

`enrichment.py` is a Python module that generates enrichment information for Next-Generation Sequencing (NGS) data. It is part of a larger NGS pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `subprocess`
- `pathlib`
- `dotenv`
- `rich`

## Functions

The script contains the following function:

- `get_enrichment()`: Generates enrichment information. It takes the sample ID, directory of the run, and panel type as arguments. It does not return anything.

## Usage

The `get_enrichment()` function is called within the larger NGS pipeline to generate enrichment information. It is not intended to be run as a standalone script.

## Output

The `get_enrichment()` function outputs an enrichment file with the enrichment information. It also logs the start of the enrichment process, whether the enrichment file exists or was created, and any errors encountered during the process.

## Example

Here is an example of how the `get_enrichment()` function might be used within the larger NGS pipeline:

```python
from src.lib.enrichment import get_enrichment

sample_id = "sample1"
datadir = "/path/to/data"
panel = "Cplus"

get_enrichment(sample_id, datadir, panel)
```

In this example, the `get_enrichment()` function is called with the sample ID, directory of the run, and panel type. The function does not return anything.