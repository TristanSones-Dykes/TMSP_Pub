# src

This directory contains helper functions and pipelines used in scripts/ to run the analysis and generate figures.

## hydrophobicity.r 

Contains functions for hydrophobicity analysis as well as an R function to clean the Phobius output.

## phobius.py

A Python script that defines `phobius`, a function for running Phobius on a proteome file (.fasta), does not require Phobius to be installed as it uses the online server.

## deepTMHMM.py

Defines `extract_deepTMHMM`, a function for extracting the results of DeepTMHMM analysis. Returning a .csv file of window types and lengths.

# psipred

Contains data and scripts for running PSIPRED on S. Cerevisiae proteins
predicted by Phobius.