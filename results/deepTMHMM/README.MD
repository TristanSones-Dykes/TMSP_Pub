This folder contains the DeepTMHMM results for S Cerevisiae, S Pombe, and Candida Albicans.

Each organism subfolder contains:
 - `deeptmhmm_results.md`: A markdown file with the counts of protein types and region types found in the entire proteome.
 - `predicted_topologies.3line`: A 3-line file containing each protein, its sequence on one line, and predicted topology on the next line.
 - `TMRs.gff3`: A GFF3 file containing the predicted transmembrane/SP regions and their start and end positions.

Data source:
 - DeepTMHMM was ran using the command: `sudo env "BIOLIB_LOG=debug" biolib run --local DTU/DeepTMHMM --fasta {input fasta directory}`, notably without
 using the gpu flag due to driver issues. If running this in future, note that it speeds up as it processes more proteins, so as long as you have enough
 cores and cpu cache it should run in a *reasonable* time.