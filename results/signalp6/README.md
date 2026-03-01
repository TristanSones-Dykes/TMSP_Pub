## signalp6

Results for running SignalP6.0 predictions.

Initial run is on S. cerevisiae proteins only.

Code used to install and run SignalP6.0, from the directory TMSP_PUB was:

```
conda create -n run-signalp6 python=3.10 numpy=1.23 pybiolib
conda activate run-signalp6 
biolib run DTU/SignalP-6 --help
mkdir results/signalp6
mkdir results/signalp6/S_Cerevisiae 
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae.fasta --output_dir results/signalp6/S_Cerevisiae --format txt --organism eukarya --mode fast
```

However, this ran via the SignalP server up to `Predicting:  12% 700/5907 [09:50<1:13:28,  1.18sequences/s]`, then gave error `biolib.biolib_errors.BioLibError: Cloud: Job exceeded max run time`

This means that, to run we will either need to break into parts...

````
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_001.fasta --output_dir results/signalp6/S_Cerevisiae_001 --format txt --organism eukarya --mode fast
# Getting FileNotFoundError: [Errno 2] No such file or directory: 'output/output.json'
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_002.fasta --output_dir results/signalp6/S_Cerevisiae_002 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_003.fasta --output_dir results/signalp6/S_Cerevisiae_003 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_004.fasta --output_dir results/signalp6/S_Cerevisiae_004 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_005.fasta --output_dir results/signalp6/S_Cerevisiae_005 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_006.fasta --output_dir results/signalp6/S_Cerevisiae_006 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_007.fasta --output_dir results/signalp6/S_Cerevisiae_007 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_008.fasta --output_dir results/signalp6/S_Cerevisiae_008 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_009.fasta --output_dir results/signalp6/S_Cerevisiae_009 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_010.fasta --output_dir results/signalp6/S_Cerevisiae_010 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_011.fasta --output_dir results/signalp6/S_Cerevisiae_011 --format txt --organism eukarya --mode fast
biolib run DTU/SignalP-6 --fastafile data/Proteins/full/S_Cerevisiae_012.fasta --output_dir results/signalp6/S_Cerevisiae_012 --format txt --organism eukarya --mode fast


````