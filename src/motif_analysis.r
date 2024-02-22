library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))
library(Biostrings)
library(ggseqlogo)
library(stringr)


# ---- Description ----
# This script contains functions for finding motifs in an RNA/DNA sequence


# function that takes a 5'UTR FASTA returns a dataframe
# with the number of motifs found in each sequence
count_motifs <- function(input_path, motif = "CNYTCNYT", isPath = TRUE) {
    # load fasta file as DNA string set
    if (isPath) {
        input_RNA <- readDNAStringSet(input_path)
    } else {
        input_RNA <- input_path
    }

    # remove truncated sequences, which we don't need.
    input_RNA <- input_RNA[width(input_RNA) == 1000]

    # extract IDs
    input_ids <- names(input_RNA)

    # find motifs
    motif_matches <- tibble(seqid = input_ids,
                            count_up_1000 = vcountPattern(pattern = DNAString(motif),
                                    subject = input_RNA,
                                    fixed = "subject"),
                            count_up_200 = vcountPattern(pattern = DNAString(motif),
                                    subject = subseq(input_RNA, start = 801L, end = 1000L),
                                    fixed = "subject"),
                            count_up_100 = vcountPattern(pattern = DNAString(motif),
                                    subject = subseq(input_RNA, start = 901L, end = 1000L),
                                    fixed = "subject"))

    # return IDs
    return(motif_matches)
}


# function that takes a 5'UTR FASTA returns a dataframe
# with the locations of motifs found in each sequence
# can use motif_lim to limit the sequences to those with
# a certain number of motifs
find_motifs <- function(input_path, motif, motif_lim = NULL, isPath = TRUE) {
    # load fasta file as DNA string set
    if (isPath) {
        input_RNA <- readDNAStringSet(input_path)
    } else {
        input_RNA <- input_path
    }

    # remove truncated sequences, which we don't need.
    input_RNA <- input_RNA[width(input_RNA) == 1000]

    # extract IDs
    input_ids <- names(input_RNA)

    # if the motif A lim is not null, calculate the motif A matches
    # and subset the input RNA by the motif A lim
    if (!is.null(motif_lim)) {
        motif_limited <- count_motifs(input_path, motif) %>%
            filter(count_up_1000 >= motif_lim) %>%
            pull(seqid)

        input_RNA <- input_RNA[motif_limited]
        input_ids <- names(input_RNA)
    }

    # find motif A
    motif_matches <- data.frame(vmatchPattern(pattern = DNAString(motif),
                                        subject = input_RNA,
                                        fixed = "subject")) %>%
                                        select(group, start, end)
    motif_matches$seqid <- input_ids[motif_matches$group]

    return(motif_matches)
}


# function that takes an IUPAC motif sequence and returns
# a PDict of all possible sequences
iupac_to_ACGT <- function(iupac_code) {
  iupac_map <- list(
    "A" = "A",
    "C" = "C",
    "G" = "G",
    "T" = "T",
    "U" = "T",
    "W" = c("A", "T"),
    "S" = c("C", "G"),
    "M" = c("A", "C"),
    "K" = c("G", "T"),
    "R" = c("A", "G"),
    "Y" = c("C", "T"),
    "B" = c("C", "G", "T"),
    "D" = c("A", "G", "T"),
    "H" = c("A", "C", "T"),
    "V" = c("A", "C", "G"),
    "N" = c("A", "C", "G", "T")
  )
  
  iupac_chars <- strsplit(iupac_code, "")[[1]]
  possible_bases <- lapply(iupac_chars, function(char) iupac_map[[char]])
  
  possible_combinations <- do.call(expand.grid, possible_bases)
  possible_sequences <- apply(possible_combinations, 1, paste, collapse = "")
  
  return(possible_sequences)
}


# function that takes a 5'UTR FASTA, finds the locations of motif A
# and searches for motif B in a range downstream of motif A
# motif A is not fixed, motif B is fixed within ambiguity codes
# returns dataframe with motif B codes and counts
count_motif_pairs <- function(input_path, motif_A = "CNYTCNYT", motif_B = "CCNNCT", window_left = 15, window_right = 4, motif_A_lim = NULL, isPath = TRUE) {
    # load fasta file and find motif A
    if (isPath) {
        input_RNA <- readDNAStringSet(input_path)
    } else {
        input_RNA <- input_path
    }
    motif_A_matches <- find_motifs(input_RNA, motif_A, motif_A_lim, isPath = FALSE)

    # subset motif A sequences by window size
    stringset <- DNAStringSet()
    for (i in seq_len(nrow(motif_A_matches))) {
        if (motif_A_matches$start[i] >= (window_left + 1)) {
            subsequence <- subseq(input_RNA[motif_A_matches$seqid[i]], start = motif_A_matches$start[i] - window_left, end = motif_A_matches$start[i] - window_right)
            stringset <- c(stringset, subsequence)
        }
    }

    # generate PDict for motif B
    motif_B_list <- iupac_to_ACGT(motif_B)

    # calculate number of matches for each motif B possibility
    motif_B_matches <- vcountPDict(PDict(motif_B_list),
                                   subject = stringset,
                                   fixed = TRUE,
                                   collapse = 1)

    # return dataframe of motif B codes and counts
    return(data.frame(
        motif_B = motif_B_list,
        count = motif_B_matches
    ))
}


# function that takes a 5'UTR FASTA, finds the locations of motif A
# and searches for motif B in a range upstream of motif A
# returns dataframe with motif A starts and motif B relative starts and B sequence
find_motif_pairs <- function(input_path, motif_A, motif_B, window_left = 15, window_right = 4, motif_A_lim = NULL, isPath = TRUE) {
    # load fasta file and find motif A
    if (isPath) {
        input_RNA <- readDNAStringSet(input_path)
    } else {
        input_RNA <- input_path
    }
    motif_A_matches <- find_motifs(input_RNA, motif_A, motif_A_lim, isPath = FALSE)
    
    # generate PDict for motif B
    motif_B_list <- iupac_to_ACGT(motif_B)
    motif_B_PDict <- PDict(motif_B_list)


    # match motif B to motif A sequences
    # add to dataframe using absolute position
    match_df <- data.frame(seqid = character(),
                           motif_A_start = integer(),
                           motif_B_start = integer(),
                           motif_B = character(),
                           stringsAsFactors = FALSE)

    for (i in seq_len(nrow(motif_A_matches))) {
        if (motif_A_matches$start[i] >= (window_left + 1)) {
            subsequence <- subseq(input_RNA[motif_A_matches$seqid[i]], start = motif_A_matches$start[i] - window_left, end = motif_A_matches$start[i] - window_right)[[1]]
            matches <- data.frame(matchPDict(motif_B_PDict, subsequence))

            for (j in seq_len(nrow(matches))) {
                # motif B is the position of A plus relative position of B
                motif_B_loc <- motif_A_matches$start[i] + (matches$start[j] - (window_left + 1))

                match_df <- rbind(match_df, data.frame(seqid = motif_A_matches$seqid[i],
                                                       motif_A_start = motif_A_matches$start[i],
                                                       motif_B_start = matches$start[j] - (window_left + 1),
                                                       motif_B = as.character(subseq(input_RNA[motif_A_matches$seqid[i]], start = motif_B_loc, end = motif_B_loc + str_length(motif_B) - 1)[[1]])))
            }
        }
    }
    
    return(match_df)
}