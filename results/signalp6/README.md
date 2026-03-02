# signalp6

Results for running SignalP6.0 predictions. Initial run is on S. cerevisiae proteins only.

# Contents

## S_Cerevisiae

SignalP 6.0i predictions on Saccharomyces cerevisiae whole proteome. Contents saved are:

- output.gff3 - predictions of entire signal peptide results in gff3 format
- region_output.gff3 - predictions of sub-regions (n-region, h-region, c-region) in gff3 format
- prediction_results.txt - table of proteome-wide results


# How SignalP was installed and run.

Code used to install and run SignalP6.0, from the directory TMSP_PUB was:

```
conda create -n run-signalp6 python=3.10 numpy=1.23 pybiolib
conda activate run-signalp6 
biolib run DTU/SignalP-6 --help
mkdir results/signalp6
mkdir results/signalp6/S_Cerevisiae 
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae.fasta --output_dir results/signalp6/S_Cerevisiae --format txt --organism eukarya --mode fast
```

However, running via biolib uses the SignalP server, has limits of < 1000 proteins, and

Installed it locally following [Installation Instructions](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md)

- Unpack the downloaded tar.gz file.
- Created the Python environment
- Open the directory containing the downloaded package, and install it by executing the following command.
- `pip install signalp-6-package/`

Then ran on a test file with 20 proteins by running:
```
signalp6 --fastafile data/Proteins/full/S_Cerevisiae_test20.fasta --output_dir results/signalp6/S_Cerevisiae_test20 --format txt --organism eukarya --mode fast --model_dir /Users/ewallac2/Downloads/signalp6_fast/signalp-6-package/models
```

This produced sensible outputs

output.gff:
```
## gff-version 3
YAL007C-t26_1	SignalP-6.0	signal_peptide	1	25	0.9997047	.	.	.
```

region_output.gff:
```
## gff-version 3
YAL007C-t26_1	SignalP-6.0	n-region	1	3	.	.	.	.
YAL007C-t26_1	SignalP-6.0	h-region	4	20	.	.	.	.
YAL007C-t26_1	SignalP-6.0	c-region	21	25	.	.	.	.
```

Ran on full set of S. cerevisiae proteins with:

```
signalp6 --fastafile data/Proteins/full/S_Cerevisiae.fasta --output_dir results/signalp6/S_Cerevisiae --format txt --organism eukarya --mode fast --model_dir /Users/ewallac2/Downloads/signalp6_fast/signalp-6-package/models
```
