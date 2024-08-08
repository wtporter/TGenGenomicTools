#' @title SNPs to Amino Acids
#'
#' @description
#' This function takes a data frame from fasta.to.snps or any data frame with a column named "SNP" and SNPS separated as "Reference_Location_SNP"
#' @param snp_db A data frame with a column named "SNP"
#' @param ref_seq A path to a reference sequence (.gb) format.
#' @param cores Option to specify number of cores to use.
#' @return a data frame with SNP columns and an AminoAcid column
#' @export snps.to.amino
#' @import genbankr
#' @import Biostrings
#' @import tidyverse
#' @import foreach
#' @import parallel
#' @import parallelly
#' @import doParallel

snps.to.amino <- function(snp_db, ref_seq, cores = parallelly::availableCores()){

  library(tidyverse)
  # Read in the reference sequence
  reference <- genbankr::readGenBank(ref_seq)
  # Save list of CDS regions
  Reference_DF <- data.frame(reference@cds)
  #Create new column of CDS nucleotide sequence
  Reference_DF$sequence <- "No Seq"
  #Fill nucleotide sequence
  for (i in 1:nrow(Reference_DF)) {
    #Extract sequence for each protein
    Reference_DF$sequence[i] <- substr(as.character(reference@sequence), Reference_DF[i,2], Reference_DF[i, 3])
  }

  #Separate SNPs
  Amino_Acid_List <- snp_db %>%
    separate(SNP, into =c("reference", "snp_position"), sep = "(?<=\\D)(?=\\d)", remove = F) %>%
    separate(snp_position, into =c("snp_position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)", remove = F)

  Amino_Acid_List$snp_position <- as.numeric(Amino_Acid_List$snp_position)

  # Register a parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  Temp <- foreach(SNP = 1:nrow(Amino_Acid_List), .combine = rbind) %dopar% {

    library(dplyr)

    REFERENCE <- Amino_Acid_List$reference[SNP]
    POSITION <- Amino_Acid_List$snp_position[SNP]
    MUTATION <- Amino_Acid_List$snp_mutation[SNP]

    Out <- data.frame()

    if(MUTATION !="_"& nchar(MUTATION) == 1){

      for (GENE in 1:nrow(Reference_DF)) {

        if( POSITION >= Reference_DF$start[GENE] & POSITION <= Reference_DF$end[GENE]) {

          SNP_in_gene <- (POSITION - Reference_DF$start[GENE])+1 #Plus 1 because if gene starts at 91, mutation is at 91 91-91=0, but the mutation is at position 1

          Reference_Seq <- Biostrings::DNAString(Reference_DF$sequence[GENE])

          Observed_Seq <- Biostrings::DNAString(Reference_DF$sequence[GENE])

          Theoretical_Ref <- Reference_Seq[SNP_in_gene]

          Observed_Seq[SNP_in_gene] <- MUTATION

          AA_Seq <- Biostrings::translate(Reference_Seq, if.fuzzy.codon = "solve")

          AA_Observed <- Biostrings::translate(Observed_Seq, if.fuzzy.codon = "solve")

          Alignment <- Biostrings::pairwiseAlignment(AA_Seq, AA_Observed)

          Mutations <- Alignment %>%
            Biostrings::mismatchTable() %>%
            mutate(AA_Change = paste(Reference_DF$gene[GENE], ":", PatternSubstring, PatternStart, SubjectSubstring, sep = ""))

          Out <- rbind(Out, data.frame("snp_position" = POSITION, "snp_mutation" = MUTATION, "Theoretical_Reference" = Theoretical_Ref, "Gene" = Reference_DF$gene[GENE], "Product" = Reference_DF$product[GENE], "AA" = ifelse(length(as.character(unique(Mutations$AA_Change))) > 0 ,
                                                                                             as.character(unique(Mutations$AA_Change)),
                                                                                             "Synonymous")))
        }
      }
    }
    Out
  }
  # Stop the parallel cluster
  stopCluster(cl)

  Temp$snp_position <- as.numeric(Temp$snp_position)

  Out <- full_join(Amino_Acid_List, Temp)

  Out$AA[is.na(Out$AA) & nchar(Out$snp_mutation) > 1] <- "Insertions Not Supported"

  Out$AA[is.na(Out$AA) & Out$snp_mutation == "_"] <- "Deletions Not Supported"

  Out$AA[is.na(Out$AA)] <- "Non-coding SNP"

  Out <- select(Out, SNP, AA, Gene, Product, Theoretical_Reference)

  Out
}

# # Example of functionality
# library(TGenGenomicTools)
# library(genbankr)
# library(Biostrings)
# library(tidyverse)
# library(foreach)
# library(parallel)
# library(parallelly)
# library(doParallel)
# snp_db <- fasta.to.snps(fasta="/scratch/tporter/RSV_20230402_SequencingMethodCompare_SQL/RSV_A_Fasta/TGen_Clinical_RSVA.fasta",
#                         ref_seq="/scratch/tporter/RSV_20230630_Phylogenetic_Analysis/",
#                         cores=parallelly::availableCores())
# snp_db <- read.csv("/scratch/tporter/RSV_20230908_HomoplasticSNPAnalysis/RSV_A_Homoplastic.csv")
# ref_seq="/scratch/tporter/RSV_20230825_ROI_Formatting/RSV_A_NC_038235.1.gb"
# snp_db <- dplyr::select(snp_db, SNP)
# snps.to.amino(snp_db, ref_seq)
