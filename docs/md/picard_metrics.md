# picard_metrics.py Documentation

## Overview

`picard_metrics.py` is a Python module that provides parsing capabilities for multiple types of Picard output. It is part of a larger Next-Generation Sequencing (NGS) pipeline.

## Dependencies

The script depends on several Python libraries, including:

- `os`
- `subprocess`
- `pathlib`
- `rich`

## Functions

The script contains the following functions:

- `get_yield()`: Generates Picard CollectQualityYieldMetrics. It takes the sample ID, directory of the run, and Picard jar file location as arguments. It returns 'success' if the Picard CollectQualityYieldMetrics file is successfully created, 'exists' if the file already exists.

- `get_hs_metrics()`: Generates Picard CollectHsMetrics. It takes the sample ID, directory of the run, reference genome, Picard specific BED file, Picard jar file location, and panel type as arguments. It returns 'success' if the Picard CollectHsMetrics file is successfully created, 'exists' if the file already exists.

- `get_align_summary()`: Generates Picard CollectAlignmentSummaryMetrics. It takes the sample ID, directory of the run, reference genome, and Picard jar file location as arguments. It returns 'success' if the Picard CollectAlignmentSummaryMetrics file is successfully created, 'exists' if the file already exists.

- `get_call_metrics()`: Generates Picard CollectVariantCallingMetrics. It takes the sample ID, directory of the run, VCF file used for variant calling, and Picard jar file location as arguments. It returns 'success' if the Picard CollectVariantCallingMetrics file is successfully created, 'exists' if the file already exists.

## Usage

The functions in `picard_metrics.py` are called within the larger NGS pipeline to generate various Picard metrics. It is not intended to be run as a standalone script.

## Output

The functions in `picard_metrics.py` output various results, including flags indicating the completion of certain functions, and files with Picard metrics.

## Example

Here is an example of how the `get_yield()` and `get_hs_metrics()` functions might be used within the larger NGS pipeline:

```python
from src.lib.picard_metrics import get_yield, get_hs_metrics

sample_id = "sample1"
datadir = "/path/to/data"
reference = "/path/to/reference"
picard = "/path/to/picard.jar"
bait_file = "/path/to/bait_file"

result_yield = get_yield(sample_id, datadir, picard)
result_hs_metrics = get_hs_metrics(sample_id, datadir, reference, bait_file, picard)
```

In this example, the `get_yield()` function is called with the sample ID, directory of the run, and Picard jar file location. The function returns 'success' if the Picard CollectQualityYieldMetrics file is successfully created, 'exists' if the file already exists. The `get_hs_metrics()` function is called with the sample ID, directory of the run, reference genome, Picard specific BED file, Picard jar file location, and panel type. The function returns 'success' if the Picard CollectHsMetrics file is successfully created, 'exists' if the file already exists.