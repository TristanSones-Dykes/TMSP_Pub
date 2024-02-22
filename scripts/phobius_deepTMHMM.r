library(here)
# grab functions from src
source(here("src", "signal_peptides.r"))
source(here("src", "hydrophobicity.r"))
source(here("src", "motif_analysis.r"))
library(mixtools)

protein_paths <- base::Map(paste, here("data", "Proteins", "full"), list.files(here("data", "Proteins", "full")), sep = "/")
species_names <- gsub(".fasta", "", list.files(here("data", "Proteins", "full")))

# run phobius
phobius_results <- lapply(protein_paths, r_phobius, fullSignal = TRUE, subset = 60)
