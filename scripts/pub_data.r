###
# pub_data.r
# The purpose of this script is to generate all data used in
# figure generation and statistical analysis for the publication.
###

library(here)
# grab functions from src
source(here("src", "hydrophobicity.r"))
source(here("src", "utils.r"))


# --- Run all sequences through phobius --- #

species_df <- read_species_table("proteome_table_fungi12.txt")
protein_paths <- here("data", "proteins", "pub", species_df$Filename)
species_names <- species_df$Nicename

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# run phobius
phobius_results <- lapply(protein_paths, run_phobius)

# filter, add window_length and species name, then combine
phobius_df <- process_phobius_results(phobius_results, species_names)

# --- Subset and create S_Cerevisiae.csv --- #

S_Cerevisiae <- phobius_df %>%
  filter(species == "Saccharomyces cerevisiae") %>%
  select(-species)

# extract Rose scale
rose <- scales[, c("aa", "Rose")]
colnames(rose) <- c("V1", "V2")

# extract Hessa scale
Hessa <- scales[, c("aa", "Hessa")]
colnames(Hessa) <- c("V1", "V2")

# read protein sequence and calculate hydropathy for KD and Rose
AA_stringset <- readAAStringSet(here("data", "proteins", "pub", "S_Cerevisiae.fasta"))
KD_df <- add_compound_hydropathy_score(S_Cerevisiae, AA_stringset, scale = KD, include_max = TRUE)
rose_df <- add_compound_hydropathy_score(S_Cerevisiae, AA_stringset, scale = rose, include_max = TRUE)
Hessa_df <- add_compound_hydropathy_score(S_Cerevisiae, AA_stringset, scale = Hessa, include_max = TRUE)

# rename and join
KD_df <- KD_df %>% rename(compound_hydropathy = "KD_hydropathy", max_hydropathy = "KD_max_hydropathy")
rose_df <- rose_df %>% rename(compound_hydropathy = "rose_hydropathy", max_hydropathy = "rose_max_hydropathy")
Hessa_df <- Hessa_df %>% rename(compound_hydropathy = "Hessa_hydropathy", max_hydropathy = "Hessa_max_hydropathy")
SC_output <- KD_df %>%
  inner_join(rose_df %>% select(seqid, rose_hydropathy, rose_max_hydropathy), by = "seqid") %>%
  inner_join(Hessa_df %>% select(seqid, Hessa_hydropathy, Hessa_max_hydropathy), by = "seqid")

# write to file
write.csv(SC_output, here("results", "figures", "SC_first_60.csv"), row.names = FALSE)

# --- Create SP and TM files for each species --- #

# write SP and TM protein ID lists for each species
for (i in seq_along(phobius_results)) {
  write_protein_lists(phobius_results[[i]], species_names[i])
}


# --- Run S_Cerevisiae through Phobius but full length results --- #

# run phobius
full_results <- run_phobius(here("data", "Proteins", "pub", "S_Cerevisiae.fasta"), fullSignal = TRUE)

# --- GO Analysis using fungidb --- #

for (i in seq_len(dim(species_df)[1])) {
  species_tax <- species_df$FungiDB_id[i]
  species_name <- species_df$Nicename[i]
  output_file <- str_split(species_df$Filename[i], "\\.")[[1]][1]

  for (type in c("TM", "SP")) {
    input_file <- here("results", "proteins", paste(species_name, paste(type, ".txt", sep = ""), sep = "_"))
    output_dir <- here("results", "GO", paste(output_file, type, sep = "_"))

    # Check if the output already exists, if not, run the GO analysis
    if (dir.exists(output_dir)) {
      cat(paste("GO analysis for", species_name, type, "proteins already exists. Skipping...\n"))
      cat("-----------------------------------\n\n")
      next
    }

    cat(paste("Running GO analysis for", species_name, type, "proteins:\n"))
    system(paste("sh", here("src", "GO_analysis.sh"), "-t", paste("'", species_tax, "'", sep = ""), "-i", paste("'", input_file, "'", sep = ""), "-o", output_dir))
    cat("-----------------------------------\n\n")
  }
}

# ran SP and TM selected gene IDs through PANTHER GO-slim for humans
# SP: Extraceullar region, TM: Membrane
human_rows <- data.frame(
  prediction = c("SP", "TM"),
  Name = c("extracellular region", "membrane"),
  p_value = c(0, 9.19e-185)
)
colnames(human_rows) <- c("prediction", "Name", "P-value")

human_rows %>%
  filter(prediction == "SP") %>%
  select(-prediction) %>%
  write.table(here("results", "GO", "human_ref_SP", "goEnrichmentResult.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
human_rows %>%
  filter(prediction == "TM") %>%
  select(-prediction) %>%
  write.table(here("results", "GO", "human_ref_TM", "goEnrichmentResult.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# --- Write S_Cerevisiae.fa for PSIPRED --- #

cerevisiae_names <- phobius_df %>%
  filter(species == "Saccharomyces cerevisiae") %>%
  pull(seqid)

library(Biostrings)

writeXStringSet(proteins[[1]][cerevisiae_names], here("src", "psipred", "S_Cerevisiae.fa"))

# --- Run human genome through Phobius --- #

full_human <- readAAStringSet(here("data", "Proteins", "protein.faa"))

# subset those that are at least 60 amino acids long
human <- full_human[width(full_human) >= 60]
human <- subseq(human, start = 1, end = 60)

# write to file
writeXStringSet(human, here("data", "Proteins", "pub", "human_ref.fasta"))

# run phobius
human_results <- run_phobius(here("data", "Proteins", "pub", "human_ref.fasta"))

# write TM and SP to file
write_protein_lists(human_results, "human")


# --- Create 60-length files for new species --- #
new_species_df <- here("data", "proteins", "full", "proteome_table.txt") %>%
  read_tsv(comment = "#") %>%
  mutate(Nicename = as_factor(Nicename))
new_protein_paths <- here("data", "proteins", "full", new_species_df$Filename)

# read in protein sequences
new_proteins <- lapply(new_protein_paths, readAAStringSet)

# remove if less than 60 amino acids
new_proteins <- lapply(new_proteins, function(x) x[width(x) >= 60])
# subset to first 60 amino acids
new_proteins <- lapply(new_proteins, function(x) subseq(x, start = 1, end = 60))

# write to pub and run phobius
subset_paths <- c()
for (i in seq_along(new_proteins)) {
  # remove the .fasta extension
  file_name <- str_split(new_species_df$Filename[i], "\\.")[[1]][1]
  # create the path
  subset_path <- here("data", "proteins", "pub", paste(file_name, "_first_60.fasta", sep = ""))
  # write the file
  writeXStringSet(new_proteins[[i]], subset_path)
  subset_paths <- c(subset_paths, subset_path)
}


# === plot full length phobius results === #
new_phobius_results <- lapply(subset_paths, run_phobius)
phobius_df <- process_phobius_results(new_phobius_results, new_species_df$Filename)
phobius_plot <- plot_helix_length_histogram(phobius_df)

# --- Plots for industry-relevant species --- #
new_species_df <- read_species_table("proteome_table_industry.txt")
new_protein_paths <- here("data", "proteins", "pub", new_species_df$Filename)

new_phobius_results <- lapply(new_protein_paths, run_phobius)
phobius_df <- process_phobius_results(new_phobius_results, new_species_df$Nicename)
phobius_plot <- plot_helix_length_histogram(phobius_df)

ggsave(here("results", "figures", "phobius_industry_species.png"),
  plot = phobius_plot,
  width = 6, height = 8, units = "in", dpi = 300
)
