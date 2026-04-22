# signalp6

Results for running SignalP6.0 predictions. Initial run is on S. cerevisiae proteins only.
We ran both "fast" and "slow_sequential" models.

# Contents

## S_Cerevisiae_fast

SignalP 6.0i fast predictions on Saccharomyces cerevisiae whole proteome. Contents saved are:

- output.gff3 - predictions of entire signal peptide results in gff3 format
- region_output.gff3 - predictions of sub-regions (n-region, h-region, c-region) in gff3 format
- prediction_results.txt - table of proteome-wide results

## S_Cerevisiae_slow_sequential

SignalP 6.0i slow sequential predictions on Saccharomyces cerevisiae whole proteome. Contents saved are:

- output.gff3 - predictions of entire signal peptide results in gff3 format
- region_output.gff3 - predictions of sub-regions (n-region, h-region, c-region) in gff3 format
- prediction_results.txt - table of proteome-wide results


# How SignalP6.0 was installed and run.

## SignalP fast

Code used to install and run SignalP6.0 fast from the repository root was:

```
conda create -n run-signalp6 python=3.10 numpy=1.23 pybiolib
conda activate run-signalp6 
biolib run DTU/SignalP-6 --help
mkdir results/signalp6
mkdir results/signalp6/S_Cerevisiae_fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae.fasta --output_dir results/signalp6/S_Cerevisiae_fast --format txt --organism eukarya --mode fast
```

However, running via biolib uses the SignalP server, has limits of < 1000 proteins, and

Installed it locally following [Installation Instructions](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md)

- Unpack the downloaded tar.gz file.
- Created the Python environment
- Open the directory containing the downloaded package, and install it by executing the following command.
- `pip install signalp-6-package/`

**Note:** replace `<path-to-signalp6-models>` with the local path to the unpacked SignalP6 model directory on your machine.

Then ran on a test file with 20 proteins by running:
```
signalp6 --fastafile data/Proteins/full/S_Cerevisiae_test20.fasta --output_dir results/signalp6/S_Cerevisiae_test20_fast --format txt --organism eukarya --mode fast --model_dir <path-to-signalp6-models>
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
signalp6 --fastafile data/Proteins/full/S_Cerevisiae.fasta --output_dir results/signalp6/S_Cerevisiae_fast --format txt --organism eukarya --mode fast --model_dir <path-to-signalp6-models>
```

and removed 1000s of `*_plot.txt` files with:
```
rm results/signalp6/S_Cerevisiae_fast/*_plot.txt
```

## SignalP slow sequential

Next downloaded slow_sequential model per SignalP instructions. Used the same signalp6 installation and con

```
conda activate run-signalp6
```

Code used to run SignalP6.0 slow_sequential on a test file from the repository root was:
```
signalp6 --fastafile data/Proteins/full/S_Cerevisiae_test20.fasta --output_dir results/signalp6/S_Cerevisiae_test20_slow --format txt --organism eukarya --mode slow-sequential --model_dir <path-to-signalp6-models>
```

Then on full proteome was:
```
signalp6 --fastafile data/Proteins/full/S_Cerevisiae.fasta --output_dir results/signalp6/S_Cerevisiae_slow --format txt --organism eukarya --mode slow-sequential --model_dir <path-to-signalp6-models>
```
