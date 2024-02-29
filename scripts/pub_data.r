###
# pub_data.r
# The purpose of this script is to generate all data used in 
# figure generation and statistical analysis for the publication.
###

library(here)
# grab functions from src
source(here("src", "signal_peptides.r"))
source(here("src", "hydrophobicity.r"))


# --- Run all sequences through phobius --- #

# read file names and remove proteome_table.txt
files <- list.files(here("data", "proteins", "pub"))
files <- files[files != "proteome_table.txt"]

# attach path to protein file names
protein_paths <- base::Map(paste, here("data", "proteins", "pub"), files, sep = "/")
species_names <- gsub(".fasta", "", files)

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# run phobius
phobius_results <- lapply(protein_paths, r_phobius)

# filter out sequences with no signal peptides
# and add species name
for (i in 1:length(phobius_results)) {
    phobius_results[[i]] <- phobius_results[[i]] %>% 
        filter(phobius_end != 0) %>%
        mutate(window_length = phobius_end - phobius_start) %>% 
        mutate(species = species_names[i])
}

# join and reset row names
phobius_df <- do.call(rbind, phobius_results)
rownames(phobius_df) <- NULL

# set species order for plots
species_order <- list("S_Cerevisiae" = 1,
              "C_Albicans" = 2,
              "N_Crassa" = 3,
              "P_Oryzae" = 4,
              "Z_Tritici" = 5,
              "A_Fumigatus" = 6,
              "S_Pombe" = 7,
              "P_Graminis" = 8,
              "U_Maydis" = 9,
              "C_Neoformans" = 10,
              "R_Delemar" = 11,
              "B_Dendrobatitis" = 12)
phobius_df$species <- factor(phobius_df$species, levels = names(species_order))


# --- Subset and create S_Cerevisiae.csv --- #

S_Cerevisiae <- phobius_df %>% 
    filter(species == "S_Cerevisiae") %>%
    select(-species)

# extract Rose scale
rose <- scales[, c("aa", "Rose")]
colnames(rose) <- c("V1", "V2")

# read protein sequence and calculate hydropathy for KD and Rose
AA_stringset <- readAAStringSet(here("data", "proteins", "pub", "S_Cerevisiae.fasta"))
KD_df <- add_compound_hydropathy_score(S_Cerevisiae, AA_stringset, useSignalP = FALSE, scale = KD, include_max = TRUE)
rose_df <- add_compound_hydropathy_score(S_Cerevisiae, AA_stringset, useSignalP = FALSE, scale = rose, include_max = TRUE)

# rename and join
KD_df <- KD_df %>% rename(compound_hydropathy = "KD_hydropathy", max_hydropathy = "KD_max_hydropathy")
rose_df <- rose_df %>% rename(compound_hydropathy = "rose_hydropathy", max_hydropathy = "rose_max_hydropathy")
SC_output <- KD_df %>% 
    inner_join(rose_df %>% select(seqid, rose_hydropathy, rose_max_hydropathy), by = "seqid")

# write to file
write.csv(SC_output, here("results", "figures", "SC_first_60.csv"), row.names = FALSE)
