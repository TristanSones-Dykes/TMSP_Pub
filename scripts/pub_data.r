###
# pub_data.r
# The purpose of this script is to generate all data used in 
# figure generation and statistical analysis for the publication.
###

library(here)
# grab functions from src
source(here("src", "hydrophobicity.r"))


# --- Run all sequences through phobius --- #

# attach path to protein file names
species_df <- here("data", "proteins", "pub", "proteome_table.txt") %>%
  read_tsv(comment = "#") %>%
  # Next line makes Nicename a factor in same order as given
  mutate(Nicename = as_factor(Nicename),
         Nicename_splitline = 
           factor(Nicename, levels = Nicename,
                  labels = str_replace(Nicename, 
                                       pattern = " ", 
                                       replacement = "\n")))
protein_paths <- here("data", "proteins", "pub", species_df$Filename)
species_names <- species_df$Nicename

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# run phobius
phobius_results <- lapply(protein_paths, run_phobius)

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

# --- Subset and create S_Cerevisiae.csv --- #

S_Cerevisiae <- phobius_df %>% 
    filter(species == "Saccharomyces cerevisiae") %>%
    select(-species)

# extract Rose scale
rose <- scales[, c("aa", "Rose")]
colnames(rose) <- c("V1", "V2")

# read protein sequence and calculate hydropathy for KD and Rose
AA_stringset <- readAAStringSet(here("data", "proteins", "pub", "S_Cerevisiae.fasta"))
KD_df <- add_compound_hydropathy_score(S_Cerevisiae, AA_stringset, scale = KD, include_max = TRUE)
rose_df <- add_compound_hydropathy_score(S_Cerevisiae, AA_stringset, scale = rose, include_max = TRUE)

# rename and join
KD_df <- KD_df %>% rename(compound_hydropathy = "KD_hydropathy", max_hydropathy = "KD_max_hydropathy")
rose_df <- rose_df %>% rename(compound_hydropathy = "rose_hydropathy", max_hydropathy = "rose_max_hydropathy")
SC_output <- KD_df %>% 
    inner_join(rose_df %>% select(seqid, rose_hydropathy, rose_max_hydropathy), by = "seqid")

# write to file
write.csv(SC_output, here("results", "figures", "SC_first_60.csv"), row.names = FALSE)

# --- Create SP and TM files for each species --- #

# write each phobius_type group of each species to a text file
for (i in 1:length(phobius_results)) {
    SP <- phobius_results[[i]] %>% 
        filter(phobius_type == "SP") %>% 
        pull(seqid)
    TM <- phobius_results[[i]] %>% 
        filter(phobius_type == "TM") %>% 
        pull(seqid)

    write.table(SP, file = paste(here("results", "proteins"), paste(species_names[i], "SP.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(TM, file = paste(here("results", "proteins"), paste(species_names[i], "TM.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}
