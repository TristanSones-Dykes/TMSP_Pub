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
upper <- 33
helix_delim <- seq(lower, upper, 10)
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
    summarise(short = sum(window_length < 13),
              long = sum(window_length >= 13)) 

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

base_ScScatMarg <- 
  ggplot(labelled_df, aes(x = window_length, y = KD_max_hydropathy, 
                                   colour = `Experimental label`, 
                                   group = `Experimental label`,
                                   size  = `Experimental label`)) +
  geom_point() + 
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
        mutate(window_length = phobius_end - phobius_start) %>% 
        mutate(species = species_df$Nicename_splitline[i])
}

# join and reset row names
phobius_df <- do.call(rbind, phobius_results)
rownames(phobius_df) <- NULL

phobius_plot <- 
  ggplot(phobius_df, aes(x = window_length, fill = phobius_type)) + 
  geom_histogram(binwidth = 1, center = 0) + 
  facet_wrap(~species, scales = "free_y", ncol = 1, 
             strip.position = "right") + 
  labs(y = "Number of proteins") + 
  scale_x_helix_length +
  scale_fill_manual("Phobius prediction", 
                    values = c("SP" = "blue", "TM" = "red")) + 
  theme(legend.position = "bottom", 
        strip.text.y.right = element_text(face = "italic", angle = 0))

# save 
ggsave(filename = here("results", "figures", "phobius_helix_length.pdf"),
       plot = phobius_plot, 
       width = 5, height = 8)

ggsave(filename = here("results", "figures", "phobius_helix_length.png"),
       plot = phobius_plot, 
       width = 5, height = 8)


# Figure 3 - DeepTMHMM

# read string of here/results/proteins/SC_deeptmhmm/predicted_topologies.3line
deeptmhmm <- read_file(here("results", "deepTMHMM", "S_Cerevisiae", "predicted_topologies.3line"))
proteins <- strsplit(deeptmhmm, ">")[[1]]

deep_df <- data.frame(seqid = character(),
                      window_type = character(),
                      window_length = numeric(),
                      protein_type = character())

for (i in 2:length(proteins)) {
    protein_str <- proteins[i]
    protein_arr <- str_split(protein_str, "\n")[[1]]

    details <- str_split(protein_arr[1], " ")[[1]]
    seqid <- details[1]
    window_type <- details[3]

    if (window_type == "GLOB") {
        next
    }

    topology_string <- protein_arr[3]

    # split topology string into groups of same characters
    topology_groups <- strsplit(topology_string, "(?<=(.))(?!\\1)", perl = TRUE)[[1]]

    if (window_type == "SP") {
        for (j in seq_along(topology_groups)) {
            if (substr(topology_groups[j], 1, 1) == "S") {
                new_row <- data.frame(seqid = seqid, window_type = window_type, window_length = str_length(topology_groups[j]), protein_type = "SP")
                deep_df <- rbind(deep_df, new_row)
            }
        }
    } else if (window_type == "TM") {
        for (j in seq_along(topology_groups)) {
			if (substr(topology_groups[j], 1, 1) == "M") {
				new_row <- data.frame(seqid = seqid, window_type = window_type, window_length = str_length(topology_groups[j]), protein_type = "TM")
                deep_df <- rbind(deep_df, new_row)
			}
		}
	} else if (window_type == "SP+TM") {
		for (j in seq_along(topology_groups)) {
			if (substr(topology_groups[j], 1, 1) == "S") {
				new_row <- data.frame(seqid = seqid, window_type = "SP", window_length = str_length(topology_groups[j]), protein_type = "SP+TM")
                deep_df <- rbind(deep_df, new_row)
			} 
            else if (substr(topology_groups[j], 1, 1) == "M") {
				new_row <- data.frame(seqid = seqid, window_type = "TM", window_length = str_length(topology_groups[j]), protein_type = "SP+TM")
                deep_df <- rbind(deep_df, new_row)
		    }
		}
	}
}

plot_df <- deep_df %>% 
    mutate(window_length = window_length - 1) %>%
    mutate(`SP_&_TM` = ifelse(protein_type == "SP+TM", "SP and TM", "SP or TM"))

plot_df %>%
    ggplot(aes(x = window_length, fill = window_type)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    labs(x = "Window Length (AA)", y = "Density", title = "DeepTMHMM window lengths") +
    scale_x_continuous(breaks = seq(0, 60, 10))

plot_df %>%
    ggplot(aes(x = window_length, fill = window_type)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100) + 
    facet_wrap(~`SP_&_TM`, scales = "free_y", ncol = 1) +
    labs(x = "Window Length (AA)", y = "Density", title = "DeepTMHMM window lengths, split into both or either SP and TM") +
    scale_x_continuous(breaks = seq(0, 60, 10))

plot_df %>%
    group_by(seqid, window_type) %>% 
    summarise(window_length = sum(window_length), `SP_&_TM` = unique(`SP_&_TM`)) %>% 
    ggplot(aes(x = window_length, fill = window_type)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100) +
    facet_wrap(~`SP_&_TM`, scales = "free_y", ncol = 1) +
    labs(x = "Window Length (AA)", y = "Density", title = "DeepTMHMM window lengths, split into both or either SP and TM, summed if multiple") +
    scale_x_continuous(breaks = seq(0, 60, 10))

plot_df %>%
    group_by(seqid, window_type) %>% 
    summarise(window_length = sum(window_length), `SP_&_TM` = unique(`SP_&_TM`)) %>% 
    ggplot(aes(x = window_length, fill = window_type)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100) +
    labs(x = "Window Length (AA)", y = "Density", title = "DeepTMHMM window lengths, summed if multiple") +
    scale_x_continuous(breaks = seq(0, 60, 10))

plot_df %>%
    group_by(seqid, window_type) %>% 
    summarise(window_length = sum(window_length), `SP_&_TM` = unique(`SP_&_TM`)) %>% 
    ggplot(aes(x = window_length, fill = window_type, colour = window_type)) +
    geom_histogram(aes(y = after_stat(density), alpha = 0.9), bins = 100, position = "identity") +
    labs(x = "Window Length (AA)", y = "Density", title = "DeepTMHMM window lengths, summed if multiple, with alpha mixing") +
    scale_x_continuous(breaks = seq(0, 60, 10))

model_df <- plot_df %>%
    group_by(seqid, window_type) %>% 
    summarise(window_length = sum(window_length), `SP_&_TM` = unique(`SP_&_TM`))

library(mixtools)
model <- normalmixEM(model_df$window_length, k = 2)
plot(model, which = 2)

labelled_df <- plot_df %>% 
    mutate(seqid = str_sub(seqid, end=-7)) %>%
    mutate(`Experimental label` = case_when(seqid %in% verified_non_srp ~ "Cleaved SP",
                           seqid %in% screened_non_srp ~ "Cleaved SP",
                           seqid %in% verified_srp ~ "Non-cleaved SP",
                           TRUE ~ "unlabelled"))

labelled_df %>% 
    group_by(seqid) %>% 
    summarise(window_length = sum(window_length), `Experimental label` = unique(`Experimental label`)) %>% 
    ggplot(aes(x = window_length, fill = `Experimental label`, colour = `Experimental label`)) +
    geom_histogram(aes(y = after_stat(density), alpha = 0.9), bins = 100, position = "identity") +
    labs(x = "Window Length (AA)", y = "Density", title = "Experimental window lengths, summed if multiple, with alpha mixing") +
    scale_x_continuous(breaks = seq(0, 60, 10))

labelled_df %>% 
    filter(`Experimental label` != "unlabelled") %>%
    group_by(seqid) %>%
    reframe(window_length = sum(window_length), `Experimental label` = unique(`Experimental label`)) %>%
    ggplot(aes(x = window_length, fill = `Experimental label`, colour = `Experimental label`)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100) +
    labs(x = "Window Length (AA)", y = "Density", title = "Experimental window lengths, summed if multiple, with alpha mixing") +
    scale_x_continuous(breaks = seq(0, 60, 10))

