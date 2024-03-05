In this description "region" refers to either a Transmembrane or Signal Peptide region.

Variables:
 - `seqid`: Sequence ID, Categorical primary ID
 - `window_length`: Length of region, Discrete
 - `window_start`: Start of region, Discrete
 - `window_end`: End of region, Discrete
 - `window_type`: Type of region, Categorical (SP/TM)
 - `window_origin`: Origin tool for region, Categorical (only Phobius currently)
 - `KD_hydropathy`: Kyte-Doolittle hydropathy score, Continuous
 - `KD_max_hydropathy`: Maximum Kyte-Doolittle hydropathy score, Continuous
 - `rose_hydropathy`: Rose hydropathy score, Continuous
 - `rose_max_hydropathy`: Maximum Rose hydropathy score, Continuous

This file is obtained by running `S_Cerevisiae.csv` through phobius and then using the found regions to calculate a local hydrophobicity score along 
the sequence (hydrophobicity window length 9 ) using the Kyte-Doolittle and Rose scales. The summary statistics are the sum and maximum value.