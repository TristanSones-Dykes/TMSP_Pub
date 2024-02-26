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
theme_set(theme_cowplot(font_size = 10))


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
upper <- 35
helix_delim <- seq(lower, upper, 10)
helix_minor <- seq(lower, upper, 5)
helix_limits <- c(lower, upper)
scale_x_helix_length <- 
  scale_x_continuous(breaks = helix_delim,
                     limits = helix_limits,
                     minor_breaks = helix_minor)

# plot Kyte-Doolittle hydrophobicity axis
rough_KD_limits = c(min(labelled_df$KD_max_hydropathy),
                    max(labelled_df$KD_max_hydropathy))
rough_KD_limits
 
scale_y_KD_hydropathy <- 
  scale_y_continuous(breaks = 0:4, 
                     limits = c(0, 4.5),
                     expand = c(0,0))

# Figure 1

# Plot of contingency table of SP/TM regions found by phobius
# and those verified experimentally

# make contingency table
contingency_table <- verified_df %>% 
    group_by(`Experimental label`) %>% 
    summarise(SP = sum(window_type == "SP"),
              TM = sum(window_type == "TM"))

contingency_table <- as.table(as.matrix(contingency_table[,2:3]))

names(dimnames(contingency_table)) <- c("Experimental label", "Phobius label")
rownames(contingency_table) <- c("SRP", "Sec63")

# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
ScHydropathy_contingency_plot <- as.grob(~vcd::mosaic(contingency_table, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))


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
  xlab("Predicted helix length (AA)") + 
  ylab("Max. hydropathy (Kyte-Doolittle)") + 
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
            rel_heights = c(0.55, 1),
            rel_widths = c(1, 0.65)
  )
ScHydropathy_scatter_marginals_plot

ggsave(filename = here("results", "figures", "ScHydropathy_scatter_marginals_KD.pdf"), 
       plot = ScHydropathy_scatter_marginals_plot, 
       width = 7, height = 6, dpi = 300)


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
protein_paths <- base::Map(paste, here("data", "proteins", "pub"), list.files(here("data", "proteins", "pub")), sep = "/")
species_names <- gsub(".fasta", "", list.files(here("data", "proteins", "pub")))

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# run phobius
phobius_results <- lapply(protein_paths, read_phobius)

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


# Figure 3 - DeepTMHMM

# read string of here/results/proteins/SC_deeptmhmm/predicted_topologies.3line
deeptmhmm <- read_file(here("results", "proteins", "SC_deeptmhmm", "predicted_topologies.3line"))
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
    summarise(window_length = sum(window_length), `Experimental label` = unique(`Experimental label`)) %>%
    ggplot(aes(x = window_length, fill = `Experimental label`, colour = `Experimental label`)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100) +
    labs(x = "Window Length (AA)", y = "Density", title = "Experimental window lengths, summed if multiple, with alpha mixing") +
    scale_x_continuous(breaks = seq(0, 60, 10))
