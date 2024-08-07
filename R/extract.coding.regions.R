#' @title Extract Coding Regions
#'
#' @description
#' This function takes in a DNAStringSet of sequences and uses a reference file to extract the coding region.
#' @param DNAStringSet A DNAStringSet object with sequences to be analyzed.
#' @param ref_seq A .gb file opened with genbankr::readGenBank to be used as the reference sequence. Check to make sure the @sequence and @transcript slots are full.
#' @param cores Option to specify number of cores to use.
#' @return A dataframe with names of the sequences and all coding regions.
#' @export extract.coding.regions
#' @import genbankr
#' @import Biostrings
#' @import tidyverse
#' @import foreach
#' @import parallel
#' @import parallelly
#' @import doParallel

extract.coding.regions <- function(DNAStringSet, ref_seq, cores = parallelly::availableCores()){

  # Assign input parameters to variables
  Genomes <- DNAStringSet
  Reference <- ref_seq

  # Register a parallel backend using doParallel
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  library(foreach)

  # Using parallel processing loop through all genomes or sequences.
  Fasta_Transcript <- foreach(GENOME = 1:length(Genomes), .combine = rbind) %dopar% {

    # Load required libraries
    library(tidyverse)
    library(foreach)

    # Extract the current DNA sequence (Fasta) from the set
    Fasta <- Genomes[GENOME]

    # Perform pairwise alignment with the reference sequence
    Alignment <- pairwiseAlignment(Reference@sequence, Fasta)

    # Extract mismatch information
    Mismatch <- Alignment %>%
      Biostrings::mismatchTable()

    # Convert reference sequences to DNAString
    Observed_Seq <- DNAString(as.character(Reference@sequence))
    Reference_Seq <- DNAString(as.character(Reference@sequence))

    # Extract transcript information from the reference sequence
    Gene_DF <- Reference@transcripts

    # Create a data frame with nucleotide mismatches
    Nuc_mismatch <- Mismatch %>%
      mutate(Nuc_Change = paste(PatternSubstring, PatternStart, SubjectSubstring, sep = ""))

    # Update "observed" sequence with nucleotide mismatches
    if(nrow(Mismatch) > 0) {
      for (i in 1:nrow(Mismatch)) {
        SNP_Loc <- Mismatch$PatternStart[i]
        SNP_Nuc <- as.character(Mismatch$SubjectSubstring[i])
        Observed_Seq[SNP_Loc] <- SNP_Nuc
      }
    }

    # Process each gene in the reference transcript table
    Out <- foreach(GENE = 1:length(Gene_DF), .combine = rbind) %do% {
      Start <- Gene_DF@ranges@start[GENE]
      End <- Gene_DF@ranges@start[GENE] + Gene_DF@ranges@width[GENE] - 1

      # Translate observed DNA sequence to amino acids
      AA_Observed <- Biostrings::translate(Observed_Seq[Start:End], if.fuzzy.codon = "solve")

      # Create a data frame with gene information
      data.frame("Name" = names(Fasta), "Gene" = Gene_DF$gene[GENE], "Product" = Gene_DF$product[GENE], "AA_Seq" = as.character(AA_Observed))
    }

    # Return the data frame for the current DNA sequence
    Out
  }

  # Stop the parallel backend
  stopCluster(cl)

  # Return the combined data frame for all DNA sequences
  return(Fasta_Transcript)
}



# # Example of functionality
# library(genbankr)
# library(Biostrings)
# library(tidyverse)
# library(foreach)
# library(parallel)
# library(parallelly)
# library(doParallel)
# #
# ref_seq="/scratch/tporter/RSV_20240215_Updated_RSV_Phylogenies/Reference/RSVA_LR699737.gb"
# DNAStringSet=readDNAStringSet("/scratch/tporter/RSV_20240215_Updated_RSV_Phylogenies/RSV_A_Typed_All_Outliersremoved2/RSV_A_Typed_All_outliersremoved2.fasta")
# DNAStringSet <- DNAStringSet[1:40]
# RSV_A_Coding <- extract.coding.regions(DNAStringSet, ref_seq, cores = parallelly::availableCores())
