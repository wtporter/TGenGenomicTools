#' @title Extract Coding Regions
#'
#' @description
#' This function takes in a DNAStringSet of sequences and uses a reference file to extract the coding region.
#' @param DNAStringSet A DNAStringSet object with sequences to be analyzed.
#' @param ref_seq A path to a reference sequence (.gb) format.
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

  # Read in the reference sequence
  Reference <- genbankr::readGenBank(ref_seq)

  Genomes <- DNAStringSet

  # Register a parallel backend using doParallel
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  Fasta_Transcript <- foreach(GENOME = 1:length(Genomes), .combine = rbind) %dopar% {
    library(tidyverse)
    library(foreach)

    Fasta <- Genomes[GENOME]

    Alignment <- pairwiseAlignment(Reference@sequence, Fasta)

    pairwiseAlignment(Reference@sequence, Fasta)

    Mismatch <- Alignment %>%
      Biostrings::mismatchTable()

    Observed_Seq <- DNAString(as.character(Reference@sequence))

    Reference_Seq <- DNAString(as.character(Reference@sequence))

    Gene_DF <- Reference@transcripts

    Nuc_mismatch <- Mismatch %>%
      mutate(Nuc_Change = paste(PatternSubstring, PatternStart, SubjectSubstring, sep = ""))

    if(nrow(Mismatch) > 0) {
      for (i in 1:nrow(Mismatch)) {

        SNP_Loc <- Mismatch$PatternStart[i]

        SNP_Nuc <- as.character(Mismatch$SubjectSubstring[i])

        #cat(paste(SNP_Loc, SNP_Nuc))

        Observed_Seq[SNP_Loc] <- SNP_Nuc
      }
    }

    Out <- foreach(GENE = 1:length(Gene_DF), .combine = rbind) %do% {

      Start <- Gene_DF@ranges@start[GENE]

      End <- Gene_DF@ranges@start[GENE]+Gene_DF@ranges@width[GENE]-1

      AA_Observed <- Biostrings::translate(Observed_Seq[Start:End], if.fuzzy.codon = "solve")

      data.frame("Name" = names(Fasta), "Gene" = Gene_DF$gene[GENE], "Product" = Gene_DF$product[GENE], "AA_Seq" = as.character(AA_Observed))
    }

    Out
  }
  # Stop the parallel backend
  stopCluster(cl)

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
#
# ref_seq="/scratch/tporter/RSV_20240215_Updated_RSV_Phylogenies/Reference/RSVA_LR699737.gb"
# DNAStringSet=readDNAStringSet("/scratch/tporter/RSV_20240215_Updated_RSV_Phylogenies/RSV_A_Typed_All_Outliersremoved2/RSV_A_Typed_All_outliersremoved2.fasta")
# DNAStringSet <- DNAStringSet[1:40]
# RSV_A_Coding <- extract.coding.regions(DNAStringSet, ref_seq)
