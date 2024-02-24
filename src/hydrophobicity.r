library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(idpr)
library(rvest)
library(reticulate)


# ---- Description ----
# This script contains functions for calculating and plotting 
# the hydrophobicity of a protein sequence(s) on a kyte-doolittle scale.
# It also includes functions for finding the hydropathy window of a protein
# using the phobius online tool (does not require it to be installed)


try(source_python(here("src", "phobius.py")), silent = TRUE)
KD <- data.frame(V1 = KDNorm$V1, V2 = round((KDNorm$V2 * 9) - 4.5, 1))

scales <- read.csv(here("data", "scales.csv"), header = TRUE)


#variables for testing
#paper_df <- read_csv(here("data", "Proteins", "tslil_paper", "tslil_paper_protein_annotations.csv"))
#accession_numbers <- paper_df %>%
#    pull(`Accession Number`) %>%
#    writeLines(here("data", "Proteins", "tslil_paper", "tslil_paper_yeast_accession.txt"))
#yeast_secretome <- readAAStringSet(here("data", "Proteins", "tslil_paper", "tslil_protein_strings.fasta"))
#example_protein <- yeast_secretome[1]
#scaledHydropathyLocal(example_protein, window = 9, scale = "Kyte-Doolittle", plotResults = TRUE)

# ---- Functions ----


# function that calculates the local hydropathy of a protein
# from idpr package, modified to use non normalised weights
hydropathy_local <- function(sequence, window = 9, scale = KD) {
    seqVector <- sequenceCheck(sequence = sequence, nonstandardResidues = c('X'), method = "stop", 
        outputType = "vector", suppressOutputMessage = TRUE, suppressAAWarning = TRUE)
    if ((window%%2) == 0) {
        stop("Window must be an odd number")
    }

    names(seqVector) <- NULL
    seqLength <- length(seqVector)
    numberResiduesAnalyzed <- seqLength - (window - 1)
    positionVector <- ((window - 1)/2 + 1):(seqLength - (window - 
        1)/2)
    centerResidueVector <- seqVector[positionVector]
    windowVector <- rep(NA, numberResiduesAnalyzed)
    scoreVector <- rep(NA, numberResiduesAnalyzed)
    for (i in seq_len(numberResiduesAnalyzed)) {
        sequenceWindow <- seqVector[i:(i + (window - 1))]
        windowVector[i] <- paste0(sequenceWindow, collapse = "")
        windowValues <- scale$V2[match(sequenceWindow, scale$V1)]
        scoreVector[i] <- sum(windowValues)/window
    }
    windowDF <- data.frame(Position = positionVector, Window = windowVector, 
        CenterResidue = centerResidueVector, WindowHydropathy = scoreVector)
        
    return(windowDF)
}
#hydropathy_local(example_protein, 9)


# function that calculates the mean hydropathy of a protein
mean_hydropathy <- function(sequence, scale = KD) {
    seqCharacterVector <- sequenceCheck(sequence = sequence, method = "stop", 
        outputType = "vector", suppressOutputMessage = TRUE)
    
    seqLength <- length(seqCharacterVector)
    scoreVector <- scale$V2[match(seqCharacterVector, scale$V1)]
    mean(scoreVector)
}
#mean_hydropathy(example_protein)


# function that uses phobius to find the hydropathy window of a protein
find_hydropathy_windows <- function(protein_AAStringSet, isString = FALSE) {

    # if it is a string, write it to a file
    if (isString) {
        writeXStringSet(protein_AAStringSet, file = "temp.fasta")
        protein_AAStringSet <- "temp.fasta"
    }
    
    # call phobius on file
    hydropathy_windows <- phobius(protein_AAStringSet)

    # if it is a string, delete the file
    if (isString) {
        file.remove("temp.fasta")
    }

    return(hydropathy_windows)
}
#find_hydropathy_windows(example_protein, isString = TRUE)
#find_hydropathy_windows(here("data", "Proteins", "tslil_paper", "ten_tslil_protein_strings.fasta"))


# function that takes in a protein AA string and returns a plot of the hydropathy window
plot_hydropathy_window <- function(protein_AAStringSet, window_size = 9) {
    # use phobius to find the hydropathy window
    hydropathy_window <- c(find_hydropathy_windows(protein_AAStringSet, isString = TRUE)[1,])

    if (hydropathy_window$start == hydropathy_window$end) {
        print("No hydropathy window found")
        return()
    }

    hydropathy_df <- hydropathy_local(protein_AAStringSet)
    max_hydropathy_row <- hydropathy_df %>%
        filter(WindowHydropathy == max(WindowHydropathy))
    max_hydropathy_row$WindowHydropathy <- round(max_hydropathy_row$WindowHydropathy, 2)

    if (dim(max_hydropathy_row)[1] > 1) {
        text_angle <- 50
    } else {
        text_angle <- 0
    }

    ggplot(hydropathy_df, aes(x = Position, y = WindowHydropathy)) + 
        geom_path() + 
        geom_hline(yintercept = mean_hydropathy(protein_AAStringSet), linetype = "dashed", colour = "grey") +
        geom_vline(xintercept = hydropathy_window$start, linetype = "dashed", colour = "red") +
        geom_vline(xintercept = hydropathy_window$end, linetype = "dashed", colour = "red") +
        geom_point(data = max_hydropathy_row, aes(x = Position, y = WindowHydropathy), colour = "black") +
        geom_text(data = max_hydropathy_row, aes(x = Position, y = WindowHydropathy, label = WindowHydropathy), vjust = -1, angle=text_angle) +
        geom_rect(aes(xmin = hydropathy_window$start, xmax = hydropathy_window$end, ymin = 0, ymax = max_hydropathy_row$WindowHydropathy[1]), alpha = 0.01, fill = "lightblue") +
        labs(x = "Position", y = "Hydropathy", title = "Hydropathy of protein", colour = "Location")
}    
#plot_hydropathy_window(example_protein, 9)


# calculates the compund hydropathy score of a protein
# by multiplying the maximum hydropathy score by the length of the window
compound_hydropathy_score <- function(protein_AAStringSet, window_size = 9) {
    hydropathy_window <- c(find_hydropathy_windows(protein_AAStringSet, isString = TRUE)[1,])
    max_hydropathy_row <- hydropathy_local(protein_AAStringSet) %>%
        filter(WindowHydropathy == max(WindowHydropathy))
    
    return(max_hydropathy_row$WindowHydropathy * ((hydropathy_window$end - hydropathy_window$start) + 1))
}
#compound_hydropathy_score(example_protein, 9)


# function that takes in a protein fasta path and returns a dataframe
# with the phobius results for each protein, if a SP or TM region
# is detected, window bounds are returned
# if fullSignal is TRUE, the full SP length is returned
r_phobius <- function(protein_AA_path, fullSignal = FALSE, subset = NULL) {
    # extract file name from path
    if (fullSignal) {
        file_name <- gsub(".fasta", "_fullSignal", basename(protein_AA_path))
    } else {
        file_name <- gsub(".fasta", "", basename(protein_AA_path))
    }
        
    # create output path
    out_path <- paste0(here("results", "phobius", file_name), ".csv")

    # check if output path exists, if it does, exit function
    if (file.exists(out_path)) {
        return(read_csv(out_path))
    }

    # check if need to subset
    if (!is.null(subset)) {
        # read in stringset and remove all strings shorter than [subset]
        AA_stringset <- readAAStringSet(protein_AA_path)
        AA_stringset <- AA_stringset[nchar(AA_stringset) >= subset]

        # select first [subset] AA of each protein
        AA_subset <- subseq(AA_stringset, start = 1, end = subset)
        writeXStringSet(AA_subset, file = here("data", "Proteins", "first_60", file_name))

        # set protein_AA_path to the new subset
        protein_AA_path <- here("data", "Proteins", "first_60", file_name)
    }

    # call phobius on file
    hydropathy_windows <- phobius(protein_AA_path, fullSignal = fullSignal) %>%
        mutate(type = case_when(
            type == "TRANSMEM" ~ "TM",
            type == "SIGNAL" ~ "SP",
            type == "NONE" ~ "OTHER"
        ))
    colnames(hydropathy_windows) <- c("seqid", "phobius_start", "phobius_end", "phobius_type")

    # write to file
    write_csv(hydropathy_windows, out_path)

    return(hydropathy_windows)
}


# a function that takes in dataframe with columns
# seqid | phobius_start | phobius_end, an AAStringSet, and 
# a hydrophobicity index and returns that dataframe with the 
# compound hydropathy score for each protein
# can optionally prioritise signalp over phobius
add_compound_hydropathy_score <- function(input_window_df, AA_stringset, scale = KD, useSignalP = FALSE, window = 9, include_max = FALSE) {
    # checking whether or not to use signalp windows
    if (useSignalP) {
        # check signalp_start and signalp_end are column names in protein_window_df
        if (!("signalp_start" %in% colnames(input_window_df))) {
            stop("signalp_start not in column names of protein_window_df")
        }
        if (!("signalp_end" %in% colnames(input_window_df))) {
            stop("signalp_end not in column names of protein_window_df")
        }

        # use signalp window if it exists, otherwise use phobius window
        protein_window_df <- input_window_df %>%
            mutate(window_start = case_when(
                !is.na(signalp_start) ~ signalp_start,
                TRUE ~ phobius_start
            )) %>%
            mutate(window_end = case_when(
                !is.na(signalp_end) ~ signalp_end,
                TRUE ~ phobius_end
            )) %>%
            mutate(window_type = case_when(
                !is.na(signalp_end) ~ signalp_type,
                TRUE ~ phobius_type
            )) %>%
            mutate(window_origin = case_when(
                !is.na(signalp_end) ~ "signalp",
                TRUE ~ "phobius"
            )) %>%
            select(-phobius_start, -phobius_end, -signalp_start, -signalp_end, -phobius_type, -signalp_type)
    } else {
        protein_window_df <- input_window_df %>%
            mutate(window_start = phobius_start) %>%
            mutate(window_end = phobius_end) %>%
            mutate(window_type = phobius_type) %>%
            mutate(window_origin = "phobius") %>%
            select(-phobius_start, -phobius_end, -signalp_start, -signalp_end, -phobius_type, -signalp_type)
    }

    # filtering out proteins with no window
    has_window <- protein_window_df %>%
        filter(window_type != "OTHER")

    # calculting compound hydropathy score for each protein
    if (include_max) {
        hydropathy_df <- data.frame(seqid = character(), compound_hydropathy = numeric(), max_hydropathy = numeric())
    } else {
        hydropathy_df <- data.frame(seqid = character(), compound_hydropathy = numeric())
    }
    for(i in seq_len(nrow(has_window))) {
        hydrophobicity_df <- hydropathy_local(AA_stringset[has_window$seqid[i]], window = window, scale = scale)
        max_hydrophobicity <- hydrophobicity_df %>%
            filter(Position >= has_window$window_start[i]) %>%
            filter(Position <= has_window$window_end[i]) %>%
            summarise(max_hydrophobicity = max(WindowHydropathy)) %>%
            pull(max_hydrophobicity)

        window_length <- has_window$window_end[i] - has_window$window_start[i] + 1
        compound_hydropathy_score <- max_hydrophobicity * window_length

        if (include_max) {
            hydropathy_df <- rbind(hydropathy_df, data.frame(seqid = has_window$seqid[i], compound_hydropathy = compound_hydropathy_score, max_hydropathy = max_hydrophobicity))
        } else {
            hydropathy_df <- rbind(hydropathy_df, data.frame(seqid = has_window$seqid[i], compound_hydropathy = compound_hydropathy_score))
        }
    }

    return(protein_window_df %>%
            full_join(hydropathy_df, by = "seqid")
    )
}
