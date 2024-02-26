library(here)
# grab functions from src
source(here("src", "signal_peptides.r"))
source(here("src", "hydrophobicity.r"))
source(here("src", "motif_analysis.r"))
library(mixtools)

###
# The purpose of this script is to extract all protein IDs from all species in the proteins file
# that are divided into two groups, cleaved and non cleaved signal peptides.
# This is using phobius labels.
#
# Also contains development versions of the figures in pub_figures.r
###


# attach path to protein file names
protein_paths <- base::Map(paste, here("data", "proteins", "pub"), list.files(here("data", "proteins", "pub")), sep = "/")
species_names <- gsub(".fasta", "", list.files(here("data", "proteins", "pub")))

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# run phobius
phobius_results <- lapply(protein_paths, r_phobius)

for (i in 1:length(phobius_results)) {
    phobius_results[[i]] <- phobius_results[[i]] %>% 
        filter(phobius_end != 0) %>%
        mutate(window_length = phobius_end - phobius_start) %>% 
        mutate(species = species_names[i])
}

# join and reset row names
phobius_df <- do.call(rbind, phobius_results)
rownames(phobius_df) <- NULL

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

plot <- ggplot(phobius_df, aes(x = window_length, colour = phobius_type)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    facet_wrap(~species, scales = "free_y", ncol = 1) + 
    labs(x = "Window Length (AA)", y = "Density") + 
    scale_x_continuous(breaks = seq(0, 35, 10))

# save 
ggsave(here("results", "figures", "phobius_window_length.jpg"), plot, width = 15, height = 50, units = "cm")

# bar plot of above categories
ggplot(phobius_df, aes(x = phobius_type, fill = phobius_type)) + 
    geom_bar(stat="count") + 
    facet_wrap(~species, scales = "free_y") + 
    labs(x = "Window Length (AA)", y = "Count")


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

# find GMM with k = 2 for S_Cerevisiae
S_Cerevisiae <- phobius_df %>% 
    filter(species == "S_Cerevisiae")

model <- normalmixEM(S_Cerevisiae$window_length, k = 2)
intersection <- uniroot(function(x) {
    dnorm(x, mean = model$mu[1], sd = model$sigma[1]) - dnorm(x, mean = model$mu[2], sd = model$sigma[2])
}, c(min(model$mu), max(model$mu)))$root

ggplot(S_Cerevisiae, aes(x = window_length)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    stat_function(fun = dnorm, args = list(mean = model$mu[1], sd = model$sigma[1]), colour = "blue") + 
    stat_function(fun = dnorm, args = list(mean = model$mu[2], sd = model$sigma[2]), colour = "green") + 
    labs(x = "Window Length (AA)", y = "Count") + 
    geom_vline(xintercept = intersection, colour = "red")

# use intersection to classify SP and TM
SP <- S_Cerevisiae %>%
    filter(window_length <= intersection) %>% 
    pull(seqid)
TM <- S_Cerevisiae %>%
    filter(window_length > intersection) %>% 
    pull(seqid)

write.table(SP, file = paste(here("results", "proteins"), paste("SC_window_length", "SP.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TM, file = paste(here("results", "proteins"), paste("SC_window_length", "TM.txt", sep = "_"), sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)


screened_non_srp <- read_tsv(here("data", "SC_screened.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)

verified_srp <- read_tsv(here("data", "SC_SRP.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)

verified_non_srp <- read_tsv(here("data", "SC_non_SRP.txt"), col_names = FALSE, show_col_types = FALSE) %>%
    pull(X1)


labelled_df <- S_Cerevisiae %>% 
    mutate(verified = case_when(seqid %in% verified_non_srp ~ "non SRP",
                           seqid %in% screened_non_srp ~ "non SRP",
                           seqid %in% verified_srp ~ "SRP",
                           TRUE ~ "unlabelled")) %>% 
    filter(verified != "unlabelled")

# make contingency table
contingency_table <- labelled_df %>% 
    group_by(verified) %>% 
    summarise(SP = sum(window_length <= intersection),
              TM = sum(window_length > intersection))

# chisq test
chisq.test(contingency_table[,2:3])

phobius_proteins_SC <- S_Cerevisiae %>%
    pull(seqid)

SC_AA <- proteins[[9]]
writeXStringSet(SC_AA[phobius_proteins_SC], file = here("results", "proteins", "SC_phobius.fasta"))

# export data to results for figures
scale_name <- "Rose"
cols <- c("aa", scale_name)
rose <- scales[, cols]
colnames(rose) <- c("V1", "V2")

protein_path <- here("data", "Proteins", "S_Cerevisiae.fasta")
AA_stringset <- readAAStringSet(protein_path)

phobius_output <- r_phobius(protein_path)
signalp_output <- signalp(protein_path)
combined_df <- signalp_output %>%
    full_join(phobius_output, by = "seqid")


# ------- #
# DEPRECATED, to be changed to allow for multiple scales at once
# ------- #
rose_df <- add_compound_hydropathy_score(combined_df, AA_stringset, useSignalP = FALSE, scale = KD, include_max = TRUE) %>% 
    drop_na(compound_hydropathy) %>% 
    mutate(window_length = window_end - window_start)

write_csv(rose_df, here("results", "figures", "S_Cerevisiae.csv"))


# -------------
# Comparing deepTMHMM to phobius
# -------------

# Comparing counts
CA_df <- phobius_df %>% 
    filter(species == "C_Albicans")

phobius_counts <- CA_df %>% 
    group_by(phobius_type) %>% 
    summarise(count = n()) %>% 
    mutate(tool = "phobius")
colnames(phobius_counts) <- c("type", "count", "tool")

# read in deepTMHMM 3line 
CA_3line <- read_delim(here("results", "deepTMHMM", "C_Albicans", "predicted_topologies.3line"), delim = "\n", col_names = FALSE)

classifications <- CA_3line[seq(1, dim(CA_3line)[1], 3),] %>% 
    mutate(div = str_locate(X1, " | ")[1]) %>% 
    mutate(X1 = str_sub(X1, div + 2, str_length(X1))) %>% 
    group_by(X1) %>% 
    summarise(count = n()) %>% 
    mutate(tool = "deepTMHMM") %>% 
    mutate(X1 = str_trim(X1))
colnames(classifications) <- c("type", "count", "tool")

# combine
combined <- rbind(phobius_counts, classifications)

no_glob <- combined %>% 
    filter(type != "GLOB")

# plot
ggplot(no_glob, aes(x = type, y = count, fill = tool)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Classification", y = "Count", title = "Comparison of phobius and deepTMHMM") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# contingency table
contingency <- combined %>% 
    filter(type != "GLOB") %>% 
    spread(tool, count) %>% 
    column_to_rownames("type")

no_glob %>% 
    group_by(tool) %>%
    summarise(total = sum(count))


# comparing total lengths
full_phobius <- r_phobius(here("data", "Proteins", "pub", "C_Albicans.fasta"), fullSignal = TRUE) %>% 
    mutate(length = phobius_end - phobius_start) %>% 
    filter(phobius_end != 0) %>% 
    mutate(SP = "Full length Phobius")

h_region_phobius <- r_phobius(here("data", "Proteins", "pub", "C_Albicans.fasta"), fullSignal = FALSE) %>% 
    mutate(length = phobius_end - phobius_start) %>% 
    filter(phobius_end != 0) %>% 
    mutate(SP = "H-region Phobius")

ggplot(full_phobius, aes(x = length, fill = phobius_type)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    labs(x = "Window Length (AA)", y = "Density") + 
    scale_x_continuous(breaks = seq(0, 35, 10))

combined_phobius <- rbind(full_phobius, h_region_phobius)

ggplot(combined_phobius, aes(x = length, fill = phobius_type)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    labs(x = "Window Length (AA)", y = "Density") + 
    scale_x_continuous(breaks = seq(0, 35, 10)) + 
    facet_wrap(~SP, scales = "free_y", ncol = 1)

# deepTMHMM
library(reticulate)
source_python(here("src", "deepTMHMM.py"))

deepTMHMM_df <- extract_deepTMHMM(here("results", "deepTMHMM", "C_Albicans", "TMRs.gff3"), 1)

deepTMHMM_df %>%
    group_by(seqid) %>% 
    summarise(count = n())

# comparing total SP length
deep_SP_window <- deepTMHMM_df %>% 
    filter(window_type == "SP") %>% 
    mutate(window_length = end - start) %>% 
    select(window_length) %>% 
    mutate(tool = "deepTMHMM")
phobius_SP_window <- full_phobius %>% 
    filter(phobius_type == "SP") %>%
    mutate(window_length = phobius_end - phobius_start) %>% 
    select(window_length) %>% 
    mutate(tool = "phobius")

deep_phobius <- rbind(deep_SP_window, phobius_SP_window)
ggplot(deep_phobius, aes(x = window_length, fill = tool)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    labs(x = "Window Length (AA)", y = "Density") + 
    scale_x_continuous(breaks = seq(0, 35, 10))

# chi-squared GOF test
counts <- deep_phobius %>% 
    group_by(tool, window_length) %>%
    summarise(count = n()) %>% 
    spread(tool, count, fill = 0)
chisq.test(counts[,2:3])


# comparing total TM length (of those with only one TM)
deep_TM_window <- deepTMHMM_df %>% 
    filter(window_type == "TM") %>% 
    mutate(window_length = end - start) %>% 
    select(window_length) %>% 
    mutate(tool = "deepTMHMM")
phobius_TM_window <- full_phobius %>%
    filter(phobius_type == "TM") %>%
    mutate(window_length = phobius_end - phobius_start) %>% 
    select(window_length) %>% 
    mutate(tool = "phobius")

deep_phobius <- rbind(deep_TM_window, phobius_TM_window)
ggplot(deep_phobius, aes(x = window_length, fill = tool)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    labs(x = "Window Length (AA)", y = "Density") + 
    scale_x_continuous(breaks = seq(0, 35, 10))

# chi-squared GOF test
counts <- deep_phobius %>% 
    group_by(tool, window_length) %>%
    summarise(count = n()) %>% 
    spread(tool, count, fill = 0)
chisq.test(counts[,2:3])

both_df <- deepTMHMM_df %>% 
    mutate(window_length = end - start) %>%
    mutate(tool = "deepTMHMM") %>%
    select(seqid, window_length, tool, protein_type, window_type) %>%
    bind_rows(full_phobius %>% 
        mutate(window_length = phobius_end - phobius_start) %>%
        mutate(window_type = phobius_type) %>%
        mutate(protein_type = phobius_type) %>%
        mutate(tool = "phobius") %>%
        select(seqid, window_length, tool, window_type, protein_type))

ggplot(both_df, aes(x = window_length, fill = tool)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    labs(x = "Window Length (AA)", y = "Density", title = "Full SP and TM length predictions by DeepTMHMM and Phobius") + 
    scale_x_continuous(breaks = seq(0, 35, 10)) + 
    facet_wrap(~window_type + protein_type, scales = "free_y", labeller = "label_both")

ggplot(both_df, aes(x = window_length, fill = tool)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    labs(x = "Window Length (AA)", y = "Density", title = "Full SP and TM length predictions by DeepTMHMM and Phobius") + 
    scale_x_continuous(breaks = seq(0, 35, 10)) + 
    facet_wrap(~window_type, scales = "free_y", ncol = 1)


# Cerevisiae
SC_TMHMM <- extract_deepTMHMM(here("results", "deepTMHMM", "S_Cerevisiae", "TMRs.gff3"), 1)
SC_phobius <- r_phobius(here("data", "Proteins", "pub", "S_Cerevisiae.fasta"), fullSignal=TRUE) %>% 
    filter(phobius_type != "OTHER") %>%
    mutate(window_length = phobius_end - phobius_start) %>%
    mutate(window_type = phobius_type) %>%
    mutate(protein_type = phobius_type) %>%
    mutate(tool = "phobius") %>% 
    select(seqid, window_length, tool, window_type, protein_type)

SC_both <- SC_TMHMM %>% 
    mutate(window_length = end - start) %>%
    mutate(tool = "deepTMHMM") %>%
    select(seqid, window_length, tool, protein_type, window_type) %>%
    rbind(SC_phobius)

ggplot(SC_both, aes(x = window_length, fill = tool)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100) +
    labs(x = "Window Length (AA)", y = "Density", title = "Full SP and TM length predictions by DeepTMHMM and Phobius") +
    scale_x_continuous(breaks = seq(0, 35, 10)) +
    facet_wrap(~window_type + protein_type, scales = "free_y", labeller = "label_both")

SC_all <- extract_deepTMHMM(here("results", "deepTMHMM", "S_Cerevisiae", "TMRs.gff3")) %>% 
    mutate(window_length = end - start) %>%
    mutate(tool = "deepTMHMM") %>%
    select(seqid, window_length, tool, protein_type, window_type) %>%
    rbind(SC_phobius)
    
SC_all %>% 
    group_by(tool, window_type) %>%
    summarise(count = n()) %>%
    spread(tool, count, fill = 0)

# S Pombe
SP_TMHMM <- extract_deepTMHMM(here("results", "deepTMHMM", "S_Pombe", "TMRs.gff3"), 1)
SP_phobius <- r_phobius(here("data", "Proteins", "pub", "S_Pombe.fasta"), fullSignal=TRUE) %>% 
    filter(phobius_type != "OTHER") %>%
    mutate(window_length = phobius_end - phobius_start) %>%
    mutate(window_type = phobius_type) %>%
    mutate(protein_type = phobius_type) %>%
    mutate(tool = "phobius") %>% 
    select(seqid, window_length, tool, window_type, protein_type)

SP_both <- SP_TMHMM %>%
    mutate(window_length = end - start) %>%
    mutate(tool = "deepTMHMM") %>%
    select(seqid, window_length, tool, protein_type, window_type) %>%
    rbind(SP_phobius)

ggplot(SP_both, aes(x = window_length, fill = tool)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100) +
    labs(x = "Window Length (AA)", y = "Density", title = "Full SP and TM length predictions by DeepTMHMM and Phobius") +
    scale_x_continuous(breaks = seq(0, 35, 10)) +
    facet_wrap(~window_type + protein_type, scales = "free_y", labeller = "label_both")
