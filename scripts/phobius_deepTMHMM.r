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

# extract deeptmhmm
deepTMHMM_paths <- base::Map(paste, here("results", "deepTMHMM"), list.files(here("results", "deepTMHMM")), "TMRs.gff3", sep = "/")

library(reticulate)
source_python(here("src", "deepTMHMM.py"))
deepTMHMM_results <- lapply(deepTMHMM_paths, extract_deepTMHMM, max_regions = 1)


# check for proteins that match and create a contingency table
# for yes, no, NA and SP, TM
# take phobius as the reference
CA <- phobius_results[[1]] %>% 
    filter(phobius_type != "OTHER") %>% 
    full_join(deepTMHMM_results[[1]], by = "seqid") %>% 
    mutate(match = ifelse(phobius_type == window_type, 1, 0)) %>% 
    select(window_type, match) %>% 
    table()

SP <- phobius_results[[2]] %>% 
    full_join(deepTMHMM_results[[2]], by = "seqid") %>% 
    mutate(match = ifelse(phobius_type == window_type, 1, 0)) %>% 
    select(window_type, match) %>% 
    table()


CA <- deepTMHMM_results[[1]] %>% 
    filter(end <= 60) %>% 
    full_join(phobius_results[[1]] %>% filter(phobius_type != "OTHER"), by = "seqid") %>% 
    mutate(match = ifelse(phobius_type == window_type, 1, 0)) %>%
    select(protein_type, match) %>%
    table()

SP <- deepTMHMM_results[[2]] %>%
    filter(end <= 60) %>%
    full_join(phobius_results[[2]] %>% filter(phobius_type != "OTHER"), by = "seqid") %>%
    mutate(match = ifelse(phobius_type == window_type, 1, 0)) %>%
    select(protein_type, match) %>%
    table()


# create contingency table
CA <- deepTMHMM_results[[1]] %>% 
    filter(end <= 60) %>%
    full_join(phobius_results[[1]] %>% filter(phobius_type != "OTHER"), by = "seqid") %>%
    group_by()