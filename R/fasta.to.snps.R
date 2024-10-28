#' @title Fasta to SNPS
#'
#' @description
#' This function takes a fasta and compares it to reference sequence and computes the SNPs within the sequence. This was built for viral genomes so may struggle with other organisms.
#' @param fasta A path to the fasta file
#' @param ref_seq A path to a reference sequence (.gb) format.
#' @param cores Option to
#' @return a data frame with all SNPs within the sequence and sequence name
#' @export fasta.to.snps
#' @import genbankr
#' @import Biostrings
#' @import tidyverse
#' @import foreach
#' @import parallel
#' @import parallelly
#' @import doParallel

# Example of functionality
# library(TGenGenomicTools)
# fasta.to.snps(fasta="/scratch/tporter/RSV_20230402_SequencingMethodCompare_SQL/RSV_A_Fasta/TGen_Clinical_RSVA.fasta",
#               ref_seq="/scratch/tporter/RSV_20230630_Phylogenetic_Analysis/RSV_A.gb",
#               cores=parallelly::availableCores())

fasta.to.snps <- function(fasta, ref_seq, cores = parallelly::availableCores()){

  library(foreach)
  #Load reference sequence
  reference <- genbankr::readGenBank(ref_seq) # Read in the reference sequence
  #Load fatsta sequence
  fasta_seqs <- Biostrings::readDNAStringSet(fasta)
  #Create parallel structure
  cl <- parallel::makeCluster(cores) #not to overload your computer
  doParallel::registerDoParallel(cl)
  #Loop through to get SNPS
  out <- foreach(FASTA = 1:length(fasta_seqs), .combine = rbind) %dopar% {
    #load tidyverse for %>%
    library(tidyverse)
    #Pairwise alignment of sequences to reference
    Alignment <- Biostrings::pairwiseAlignment(pattern = reference@sequence, subject = Biostrings::subseq(fasta_seqs)[FASTA])
    #Get SNPS from alignment
    Alignment %>%
      Biostrings::mismatchTable() %>%
      filter(SubjectSubstring != "N") %>%
      mutate(SNP = paste(PatternSubstring, PatternStart, SubjectSubstring, sep = "_")) %>%
      mutate(Fasta = names(fasta_seqs)[FASTA]) %>%
      select(Fasta, SNP)
  }
  # return out
  return(out)
}



