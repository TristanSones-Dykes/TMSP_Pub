# TMSP_Pub

Analysis of protein secretion in fungi, focusing on hydrophobic regions of N-terminal transmembrane regions and signal peptides.

This repository accompanies the in-progress manuscript:

> Protein secretion routes in fungi are mostly determined by the length of the hydrophobic helix in the signal peptide.
>
> Tristan Sones-Dykes, Atlanta Cook, Edward Wallace. 

The analysis here was performed by:

- Tristan Sones-Dykes, @TristanSones-Dykes, tsd5@st-andrews.ac.uk
- Edward Wallace, @ewallace, Edward.Wallace@ed.ac.uk

# Contents

See README files in individual directories for details.

## data

Input data used for analysis.

- `Proteins` - Protein sequences in .fasta format.
- Lists of S. cerevisiae genes that use known secretion pathways in .txt format
- Hydrophobicity scales for amino acids in .csv format.

## results

All results of data analysis.

- Signal peptide and transmembrane region predictions by Phobius, in a variety of data formats
- Signal peptide and transmembrane region predictions by DeepTMHMM, in a variety of data formats
- Reports of data analysis in .html format
- Figures produced by data analysis


## scripts

Scripts for producing figures and reports (.Rmd, .R)

## src

Source scripts for running protein-level analysis (.py, .R)
