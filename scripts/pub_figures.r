###
# pub_figures.r
# The purpose of this script is to generate all figures for the writeup
# It uses post-processed data in results so that it can be run without
# signalP and phobius related dependencies.
###

library(here)
library(tidyverse)
library(cowplot)
library(ggExtra)
library(ggplotify)
library(ggthemes)
library(vcd)
library(Biostrings)
theme_set(theme_cowplot(font_size = 10) + 
            theme(strip.background = element_blank(),
                  plot.margin = unit(c(0,0,0,0), units = "inches")))


# Load data with S. cerevisiae hydropathy calculation

hydropathy_df <- read_csv(here("results", "figures", "SC_first_60.csv"))

screened_non_srp <- read_lines(here("data", "SC_screened.txt"))
verified_srp     <- read_lines(here("data", "SC_SRP.txt"))
verified_non_srp <- read_lines(here("data", "SC_non_SRP.txt"))

labelled_df <- hydropathy_df %>% 
  mutate(`Experimental label` = 
           case_when(seqid %in% verified_non_srp ~ "Sec63-dependent",
                     seqid %in% screened_non_srp ~ "Sec63-dependent",
                     seqid %in% verified_srp ~ "SRP-dependent",
                     TRUE ~ "Unverified"))

verified_df <- labelled_df %>% 
  filter(`Experimental label` != "Unverified")

# plot helix length axis
lower <- 5 
upper <- 34
helix_delim <- seq(lower, upper, 5)
helix_minor <- seq(lower, upper, 5)
helix_limits <- c(lower, upper)
scale_x_helix_length <- 
  scale_x_continuous("Predicted helix length (AA)",
                     breaks = helix_delim,
                     limits = helix_limits,
                     minor_breaks = helix_minor,
                     expand = expansion(mult = 0, add = 0.6))

# plot Kyte-Doolittle hydrophobicity axis
rough_KD_limits = c(min(labelled_df$KD_max_hydropathy),
                    max(labelled_df$KD_max_hydropathy))
rough_KD_limits
 
scale_y_KD_hydropathy <- 
  scale_y_continuous("Max. hydropathy (Kyte-Doolittle)",
                     breaks = 0:4, 
                     limits = c(0, 4.5),
                     expand = c(0,0))

# Figure 1

# Plot of contingency table of helix length predicted by phobius
# and those verified experimentally

# make contingency table for display
contingency_df_length <- verified_df %>% 
    mutate(`Experimental label` = 
           factor(`Experimental label`,
                  levels = c("Sec63-dependent", "SRP-dependent"))
    ) %>%
    group_by(`Experimental label`) %>% 
    summarise(short = sum(window_length < 14),
              long = sum(window_length >= 14)) 

contingency_table_length <- as.table(as.matrix(contingency_df_length[,2:3]))

# make table pretty for display
names(dimnames(contingency_table_length)) <- c("Experimental label", "Helix length")
rownames(contingency_table_length) <- c("Sec63", "SRP")
# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
ScHydropathy_contingency_plot_length <- as.grob(~vcd::mosaic(contingency_table_length, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))


# Contingency table of SP/TM regions predicted by phobius
# and those verified experimentally

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_length <- chisq.test(contingency_table_length)
p_value_length <- test_length$p.value
p_value_length

# make contingency table for display
contingency_df_label <- verified_df %>% 
  mutate(`Experimental label` = 
           factor(`Experimental label`,
                  levels = c("Sec63-dependent", "SRP-dependent"))
  ) %>%
  group_by(`Experimental label`) %>% 
  summarise(SP = sum(window_type == "SP"),
            TM = sum(window_type == "TM")) 

contingency_table_label <- as.table(as.matrix(contingency_df_label[,2:3]))

# make table pretty for display
names(dimnames(contingency_table_label)) <- c("Experimental label", "Phobius label")
rownames(contingency_table_label) <- c("Sec63", "SRP")

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_label <- chisq.test(contingency_table_label)
p_value_label <- test_label$p.value



# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
ScHydropathy_contingency_plot_label <- as.grob(~vcd::mosaic(contingency_table_label, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))



# Plot of contingency table of compound hydropathy
# and those verified experimentally

# make contingency table for display
contingency_df_compound <- verified_df %>% 
  mutate(`Experimental label` = 
           factor(`Experimental label`,
                  levels = c("Sec63-dependent", "SRP-dependent"))
  ) %>%
  mutate(compound = window_length * KD_max_hydropathy) %>%
  group_by(`Experimental label`) %>% 
  summarise(hi = sum(compound < 40),
            lo = sum(compound >= 40)) 

contingency_table_compound <- as.table(as.matrix(contingency_df_compound[,2:3]))

# make table pretty for display
names(dimnames(contingency_table_compound)) <- c("Experimental label", "Compound hydropathy")
rownames(contingency_table_compound) <- c("Sec63", "SRP")

# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
ScHydropathy_contingency_plot_compound <- as.grob(~vcd::mosaic(contingency_table_compound, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_compound <- chisq.test(contingency_table_compound)
p_value_compound <- test_compound$p.value
p_value_compound



# Scatter plot of all proteins with SP/TM regions found by phobius
# Highlighting verified proteins

breaks_explabel <- c("Sec63-dependent", 
                     "SRP-dependent",
                     "Unverified")
colour_explabel <- c("Sec63-dependent" = "blue", 
                     "SRP-dependent" = "red",
                     "Unverified" = "purple")
size_explabel <- c("Sec63-dependent" = 1.5, 
                   "SRP-dependent" = 1.5,
                   "Unverified" = 0.5)

compound_hydropathy_40linedf <-
  tibble(window_length = seq(4,34, 0.2),
         KD_max_hydropathy = 40 / window_length)

base_ScScatMarg <- 
  ggplot(labelled_df, aes(x = window_length, 
                          y = KD_max_hydropathy)) +
  geom_point(aes(colour = `Experimental label`, 
                 group  = `Experimental label`,
                 size   = `Experimental label`)) + 
  geom_vline(xintercept = 13.5, linetype = "dashed") + 
  geom_line(data = compound_hydropathy_40linedf,
            linetype = "dotted") + 
  # ggtitle("Phobius detected SP/TM regions") + 
  scale_x_helix_length +
  scale_y_KD_hydropathy + 
  scale_colour_manual(breaks = breaks_explabel, values = colour_explabel) +
  scale_size_manual(breaks = breaks_explabel,values = size_explabel) +
  theme(legend.box.background = element_rect(colour = "grey60"),
        legend.box.margin = margin(2, 2, 2, 2, unit = "pt"),
        legend.justification = c(0.5, 0.5),
        plot.margin = unit(c(.1,.1,.1,.1), "mm"))

side_ScScatMarg <- 
  ggplot(labelled_df  %>% 
           dplyr::mutate(`Experimental label` = 
                           factor(`Experimental label`,
                                  levels = rev(breaks_explabel))), 
         aes(y = KD_max_hydropathy,
             fill = `Experimental label`,
             group = `Experimental label`)) +
  geom_histogram(binwidth = 0.2) +
  scale_y_KD_hydropathy + 
  scale_fill_manual(breaks = breaks_explabel, values = colour_explabel) +
  facet_grid(.~ `Experimental label`, scales = "free_x") +
  labs(x = "Number of proteins") + 
  theme(legend.position = "none", 
        strip.background = element_blank(),
        strip.text = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "mm"))

top_ScScatMarg <- ggplot(labelled_df %>% 
                   dplyr::mutate(`Experimental label` = 
                                   factor(`Experimental label`,
                                          levels = breaks_explabel)), 
                 aes(x = window_length, 
                     fill = `Experimental label`, 
                     group = `Experimental label`)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 13.5, linetype = "dashed") + 
  scale_x_helix_length + 
  scale_fill_manual(values = colour_explabel) +
  facet_grid(`Experimental label` ~., scales = "free_y") +
  labs(y = "Number of proteins") + 
  theme(legend.position = "none", 
        strip.background = element_blank(),
        strip.text = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "mm"))

ScHydropathy_scatter_marginals_plot <- 
  plot_grid(top_ScScatMarg, 
            get_legend(base_ScScatMarg),
            base_ScScatMarg + 
              theme(legend.position = "none"), 
            side_ScScatMarg,
            ncol = 2,
            align = "hv",
            axis = "bl",
            rel_heights = c(0.6, 1),
            rel_widths = c(1, 0.75)
  )
ScHydropathy_scatter_marginals_plot

# save plot
ggsave(filename = here("results", "figures", "ScHydropathy_scatter_marginals_KD.pdf"), 
       plot = ScHydropathy_scatter_marginals_plot, 
       width = 6.5, height = 5.5, dpi = 300)

# save plot as .png for google docs input
ggsave(filename = here("results", "figures", "ScHydropathy_scatter_marginals_KD.png"), 
       plot = ScHydropathy_scatter_marginals_plot, 
       width = 6.5, height = 5.5, dpi = 300)

# combine figures into grid - not yet done
# Fig_1 <- plot_grid(Fig_1A, Fig_1B, Fig_1C, Fig_1D, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)

# save figure
# ggsave(here("results", "figures", "Fig_1.jpg"), Fig_1, width = 15, height = 10, dpi = 300)


# Figure 2 - histograms of window lengths for each species

read_phobius <- function(protein_AA_path) {
    # extract file name from path, replace .fasta with _out
    file_name <- gsub(".fasta", "", basename(protein_AA_path))

    # create output path
    out_path <- paste0(here("results", "phobius", file_name), ".csv")

    # check if output path exists, if it does, exit function
    if (file.exists(out_path)) {
        return(read_csv(out_path))
    } else {
        cat("No output file found")
    }
}

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

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# run phobius
phobius_results <- lapply(protein_paths, read_phobius)

for (i in 1:length(phobius_results)) {
    phobius_results[[i]] <- phobius_results[[i]] %>% 
        filter(phobius_end != 0) %>%
        mutate(window_length = phobius_end - phobius_start + 1) %>% 
        mutate(species = species_df$Nicename_splitline[i])
}

# join and reset row names
phobius_df <- do.call(rbind, phobius_results)
rownames(phobius_df) <- NULL

# get GO analysis values
species_file_names <- lapply(species_df$Filename, function(x) gsub(".fasta", "", x))

GO_df <- data.frame(species = character(),
                    prediction = character(),
                    GO_term = character(),
                    p_value = numeric())
for (species_file in species_file_names) {
    for (prediction in c("SP", "TM")) {
        prediction_file <- paste(species_file, prediction, sep = "_")
        path <- here("results", "GO", prediction_file, "goEnrichmentResult.tsv")
        
        # get lowest p-value GO term
        # after removing meaningless "cellular component" 
        go_df <- read_tsv(path)
        go_df <- go_df %>% 
            filter(Name != "cellular component") %>%
            filter(`P-value` == min(`P-value`))

        GO_df <- rbind(GO_df, 
                               data.frame(species = species_file,
                                          prediction = prediction,
                                          GO_term = go_df$Name,
                                          p_value = go_df$`P-value`))
    }
}

# heights of histogram modes for GO label positions
sub_figure_heights <- phobius_df %>% 
    filter(window_length == 12) %>%
    group_by(species) %>% 
    summarise(height = n())

# create a dataframe with species, prediction, GO term and p-value
GO_summary_df <- species_df %>% 
    mutate(Filename = gsub(".fasta", "", Filename)) %>% 
    left_join(GO_df, by = c("Filename" = "species")) %>% 
    select(species = Nicename_splitline, prediction, GO_term, p_value) %>% 
    mutate(pred_substr = paste(prediction, ": ", GO_term, sep = "")) %>%
    group_by(species) %>%
    summarise(GO_term = paste(pred_substr, collapse = "\n")) %>% 
    # shorten longest GO term for display
    mutate(GO_term = stringr::str_remove(GO_term, pattern = "external ")) %>%
    left_join(sub_figure_heights, by = c("species" = "species"))

phobius_plot <- 
  ggplot(phobius_df, aes(x = window_length, fill = phobius_type)) + 
  geom_histogram(binwidth = 1, center = 0) + 
  geom_vline(xintercept = 13.5, linetype = "dashed") +
  geom_label(data = GO_summary_df, 
             aes(x = 28, y = height %/% 1.5, label = GO_term), 
             size = 2, inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(~species, scales = "free_y", ncol = 1, 
             strip.position = "left") + 
  scale_y_continuous("Number of proteins", position = "right") + 
  scale_x_helix_length +
  scale_fill_manual("Phobius prediction", 
                    values = c("SP" = "blue", "TM" = "red")) + 
  theme(legend.position = "bottom", 
        strip.text.y.left = element_text(face = "italic", angle = 0),
        strip.placement = "outside")

# Make fungal species tree / cladogram to inform phobius plot
library(treeio)
library(ggtree)

# first define the tree in newick format, read in to treeio format
fungal12tree_data <- 
  "((((((((Sc:1,Ca:1):1,((Nc:1,Mg:1):1,(Zt:1,Af:1):1):1),Sp:1),(Pg:1,(Um:1,Cn:1):1):1):1):1,Rd:1):1,Bd:1):1);" %>%
  textConnection() %>%
  read.newick()

# Plot the tree using ggtree.
# ladderize = FALSE preserves input tip order.
fungal12tree_plot <- 
  ggtree(fungal12tree_data,
       ladderize = FALSE) +
  # geom_tiplab here would print the tip labels, useful for checking.
  # geom_tiplab() +
  scale_y_reverse()

fungal12tree_plot

# make composite plot, including moving the y-axis title 
phobius_composite_plot <- 
  plot_grid(
    fungal12tree_plot +
      theme(plot.margin = margin(t = 0, r = 0, b = 0.55, l = 0, unit = "in")),
    phobius_plot,
    nrow = 1,
    rel_widths = c(0.2,1)) + 
  theme(plot.background = element_rect(fill = "white", colour = NA))

# save 
ggsave(filename = here("results", "figures", "phobius_helix_length.pdf"),
       plot = phobius_composite_plot, 
       width = 6, height = 8)

ggsave(filename = here("results", "figures", "phobius_helix_length.png"),
       plot = phobius_composite_plot, 
       width = 6, height = 8)


# Figure 3 - DeepTMHMM

# read string of here/results/proteins/SC_deeptmhmm/predicted_topologies.3line
deeptmhmm_3line <- 
  here("results", "deepTMHMM", "S_Cerevisiae", "predicted_topologies.3line") %>%
  read_file() %>%
  # remove irrelevant transcript id
  str_remove_all("-t26_1") %>%
  strsplit(split = ">") %>%
  .[[1]]

deeptmhmm_df <- data.frame(seqid = character(),
                           DeepTMHMM_type = character(),
                           DeepTMHMM_length = numeric())

# extract SP and TM regions, within first 60 AA
for (i in 2:length(deeptmhmm_3line)) {
    # split 3line into array of strings
    protein_str <- deeptmhmm_3line[i]
    protein_arr <- str_split(protein_str, "\n")[[1]]

    # extract protein details
    details <- str_split(protein_arr[1], " ")[[1]]
    seqid <- details[1]
    DeepTMHMM_type <- details[3]

    # skip proteins that are not TM or SP
    if (DeepTMHMM_type == "GLOB") {
        next
    }

    # split topology string into groups of same characters
    topology_string <- protein_arr[3]
    topology_groups <- strsplit(topology_string, "(?<=(.))(?!\\1)", perl = TRUE)[[1]]

    # find end positions and filter for those <= 60 AA
    end_positions <- cumsum(lapply(topology_groups, str_length))
    topology_groups <- topology_groups[end_positions <= 60]

    # populate dataframes with regions
    SP_seqid_list <- c()
    for (j in seq_along(topology_groups)) {
        if (substr(topology_groups[j], 1, 1) == "S") {
            new_row <- data.frame(seqid = seqid, DeepTMHMM_type = "SP", DeepTMHMM_length = str_length(topology_groups[j]))
            deeptmhmm_df <- rbind(deeptmhmm_df, new_row)
            break
        } 
        else if (substr(topology_groups[j], 1, 1) == "M") {
            new_row <- data.frame(seqid = seqid, DeepTMHMM_type = "TM", DeepTMHMM_length = str_length(topology_groups[j]))
            deeptmhmm_df <- rbind(deeptmhmm_df, new_row)
            break
        }
    }
}


# --- Comparing to phobius ---#

# read full length phobius results
full_phobius <- read_csv(here("results", "phobius", "S_Cerevisiae_fullSignal.csv")) %>% 
    filter(phobius_end != 0) %>%
    mutate(phobius_length = phobius_end - phobius_start + 1) %>% 
    select(seqid, Phobius_type = phobius_type, Phobius_length = phobius_length) %>% 
    mutate(method = "Phobius")

labelled_phobius <- full_phobius %>% 
  mutate(`Experimental label` = 
           case_when(seqid %in% verified_non_srp ~ "Sec63-dependent",
                     seqid %in% screened_non_srp ~ "Sec63-dependent",
                     seqid %in% verified_srp ~ "SRP-dependent",
                     TRUE ~ "Unverified"))
verified_phobius <- labelled_phobius %>%
    filter(`Experimental label` != "Unverified")

# combined table of labelled proteins by classification method
combined_labelled <- full_join(labelled_phobius, deeptmhmm_df, by = "seqid")

# Compare DeepTMHMM predictions with verified proteins
# make contingency table for display
contingency_df_deeplabel <- 
  combined_labelled %>%
  filter(`Experimental label` != "Unverified") %>% 
  mutate(`Experimental label` = 
           factor(`Experimental label`,
                  levels = c("Sec63-dependent", "SRP-dependent"))
  ) %>%
  group_by(`Experimental label`) %>% 
  summarise(SP = sum(DeepTMHMM_type == "SP", na.rm = TRUE),
            TM = sum(DeepTMHMM_type == "TM", na.rm = TRUE)) 

contingency_table_deeplabel <- as.table(as.matrix(contingency_df_deeplabel[,2:3]))

# make table pretty for display
names(dimnames(contingency_table_deeplabel)) <- c("Experimental label", "DeepTMHMM label")
rownames(contingency_table_deeplabel) <- c("Sec63", "SRP")

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_deeplabel <- chisq.test(contingency_table_deeplabel)
p_value_deeplabel <- test_deeplabel$p.value
p_value_deeplabel

knitr::kable(contingency_table_deeplabel, caption = "Contingency table of SP/TM regions predicted by DeepTMHMM and verified translocon", format = "simple")



# contingency table by method
contingency_df_deepphob <- 
  combined_labelled %>%
  group_by(DeepTMHMM_type, Phobius_type) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  pivot_wider(names_from = Phobius_type, values_from = count, values_fill = 0)

contingency_table_deepphob <- as.table(as.matrix(contingency_df_deepphob[1:2,2:3]))

# make table pretty for display
names(dimnames(contingency_table_deepphob)) <- c("DeepTMHMM", "Phobius")
rownames(contingency_table_deepphob) <- c("SP", "TM")

# test if Phobius and DeepTMHMM agree, where they both make predictions
# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_deepphob <- chisq.test(contingency_table_deepphob)
p_value_deepphob <- test_length$p.value
p_value_deepphob

knitr::kable(contingency_table_deepphob, caption = "Contingency table of SP/TM regions predicted by both DeepTMHMM and Phobius", format = "simple")

# scatter of predicted lengths by method, coloured by label match

deepphob_lengthcor_df <- 
  combined_labelled %>%
  dplyr::rename(Phobius = Phobius_type, DeepTMHMM = DeepTMHMM_type) %>%
  group_by(Phobius, DeepTMHMM) %>%
  summarise(cor_length = cor(Phobius_length, DeepTMHMM_length, 
                             use = "pairwise.complete.obs"),
            count = n(),
            .groups = "drop") %>%
  mutate(cor_label = paste0("R = ", round(cor_length, digits = 2)),
         count_label = paste0("n = ", count),
         both_label = paste0(count_label, "\n", cor_label))

deepphob_match_df <- 
  combined_labelled %>%
  dplyr::rename(Phobius = Phobius_type, DeepTMHMM = DeepTMHMM_type) %>%
  group_by(Phobius, Phobius_length,  DeepTMHMM, DeepTMHMM_length) %>%
  summarise(count = n(), .groups = "drop")

deepphob_match_plot <- 
  ggplot(data = deepphob_match_df  %>%
           drop_na() %>%
           dplyr::mutate(DeepTMHMM = factor(DeepTMHMM, 
                                            levels = c("TM", "SP"))),
         aes(x = Phobius_length, y = DeepTMHMM_length)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") + 
  geom_point(aes(size = count), colour = "forestgreen") +
  geom_text(data = deepphob_lengthcor_df %>%
              drop_na() %>%
              dplyr::mutate(DeepTMHMM = factor(DeepTMHMM, 
                                               levels = c("TM", "SP"))),
            aes(label = both_label),
            x = 40, y = 60, hjust = 0, vjust = 1, size = 3,
            inherit.aes = FALSE) +
  facet_grid(DeepTMHMM ~ Phobius, labeller = label_both) +
  theme(panel.border = element_rect(fill = NA, colour = "grey90")) + 
  tune::coord_obs_pred() + 
  labs(x = "Phobius predicted length",
       y = "DeepTMHMM predicted length")

deepphob_match_plot

# save plot
ggsave(filename = here("results", "figures", "Phobius_DeepTMHMM_length_match.pdf"), 
       plot = deepphob_match_plot, 
       width = 6.5, height = 5.5, dpi = 300)

ggsave(filename = here("results", "figures", "Phobius_DeepTMHMM_length_match.png"), 
       plot = deepphob_match_plot, 
       width = 6.5, height = 5.5, dpi = 300)
 