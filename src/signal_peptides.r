library(Biostrings)
library(rtracklayer)

# ---- Description ----
# This script contains functions for finding 
# signal peptides in a protein sequence using signalp
# command line tool (requires signalp6 to be installed)
# instructions: https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md

signalp <- function(input_path, organism = "euk", mode = "fast", format = "none") {
    # extract file name from path, replace .fasta with _out
    file_name <- gsub(".fasta", "", basename(input_path))

    # create output path
    out_path <- here("results", file_name)

    # check if output path exists, if it does, exit function
    if (!file.exists(out_path)) {
        # call signalp on file
        system(paste("signalp6 -ff", input_path, "-org", organism, "-od", out_path, "-fmt", format, "--mode", mode))
    }

    # read signalp output
    signalp_output <- read.table(paste0(out_path, "/prediction_results.txt"), header = FALSE, sep = "\t")
    colnames(signalp_output) <- c("seqid", "signalp_type", "p_other", "p_SP", "CS_pos")
    signalp_output <- signalp_output %>% select(seqid, signalp_type, p_other, p_SP)

    region_output <- readGFF(paste0(out_path, "/region_output.gff3")) %>%
        filter(type == "h-region") %>%
        dplyr::rename("signalp_start" = "start", "signalp_end" = "end") %>%
        select(seqid, signalp_start, signalp_end)

    # merge signalp and region output
    final_output <- signalp_output %>% 
        left_join(region_output, by = "seqid")

    return(final_output)
}
