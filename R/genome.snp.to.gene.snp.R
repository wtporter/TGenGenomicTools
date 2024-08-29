#' @title SNPs to Amino Acids
#'
#' @description
#' This function takes a data frame from fasta.to.snps or any data frame with a column named "SNP" and SNPS separated as "Reference_Location_SNP" and uses a Genbank file to convert genome SNPs to appropriate gene snps. The function will handle
#' @param snp_db A data frame with a column named "SNP"
#' @param ref_seq A path to a reference sequence (.gb) format.
#' @param cores Option to specify number of cores to use.
#' @return a data frame with SNP columns and associated gene SNP
#' @export genome.snp.to.gene.snp
#' @import genbankr
#' @import Biostrings
#' @import tidyverse
#' @import foreach
#' @import parallel
#' @import parallelly
#' @import doParallel

genome.snp.to.gene.snp <- function(snp_db, ref_seq, cores = parallelly::availableCores()){

  library(tidyverse)
  # Read in the reference sequence
  reference <- genbankr::readGenBank(ref_seq)
  # Save list of CDS regions
  Reference_DF <- left_join(data.frame(reference@genes), data.frame(reference@cds) %>%
                              select(locus_tag, product, translation))

  #Create new column of CDS nucleotide sequence
  Reference_DF$sequence <- "No Seq"
  #Fill nucleotide sequence
  for (i in 1:nrow(Reference_DF)) {
    #Extract sequence for each protein
    Reference_DF$sequence[i] <- substr(as.character(reference@sequence), Reference_DF[i,2], Reference_DF[i, 3])
  }

  Reference_DF <- Reference_DF %>%
    mutate(gene=ifelse(is.na(gene), locus_tag, gene))

  #Separate SNPs
  SNP_List <- snp_db %>%
    separate(SNP, into =c("reference", "snp_position"), sep = "(?<=\\D)(?=\\d)", remove = F) %>%
    separate(snp_position, into =c("snp_position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)", remove = F)

  SNP_List$snp_position <- as.numeric(SNP_List$snp_position)

  # Register a parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  Temp <- foreach(SNP = 1:nrow(SNP_List), .combine = rbind) %dopar% {

    library(dplyr)

    REFERENCE <- SNP_List$reference[SNP]
    POSITION <- SNP_List$snp_position[SNP]
    MUTATION <- SNP_List$snp_mutation[SNP]

    GENOME_SNP <- paste0(REFERENCE, POSITION,MUTATION)

    Out <- data.frame()

    #if(MUTATION != "_" & nchar(MUTATION) == 1){

      for (GENE in 1:nrow(Reference_DF)) {

        if( POSITION >= Reference_DF$start[GENE] & POSITION <= Reference_DF$end[GENE]) {

            SNP_in_gene <- (POSITION - Reference_DF$start[GENE]) + 1 #Plus 1 because if gene starts at 91, mutation is at 91 91-91=0, but the mutation is at position 1

            Reference_Seq <- Biostrings::DNAString(Reference_DF$sequence[GENE])

            if(MUTATION != "_"){

              Observed_Seq <- Biostrings::DNAString(MUTATION)

            } else {

              Observed_Seq <- "_"
            }

            Theoretical_Ref <- Reference_Seq[SNP_in_gene]

            #Observed_Seq[SNP_in_gene] <- MUTATION

            #Observed_SNP <- Observed_Seq[SNP_in_gene]

            ## Add functionality to allow SNPs on second strand.
            if(Reference_DF$strand[GENE] == "-" ){

              SNP_in_gene <- (Reference_DF$end[GENE] - POSITION) + 1

              Reference_Seq <- Biostrings::reverseComplement(Reference_Seq)

              if(MUTATION != "_"){

                Observed_Seq <- Biostrings::reverseComplement(Observed_Seq)

              } else {

                Observed_Seq <- "_"
              }

              Theoretical_Ref <- Biostrings::reverseComplement(Theoretical_Ref)

              MUTATION <- Observed_Seq

            }

              Out <- rbind(Out, data.frame("SNP" = GENOME_SNP,
                                         "snp_position_genome" = POSITION,
                                         "snp_position_gene" = SNP_in_gene,
                                         "snp_mutation" = MUTATION,
                                         "Theoretical_Reference" = Theoretical_Ref,
                                         "Gene" = Reference_DF$gene[GENE],
                                         "SNP_Gene" = paste0(Theoretical_Ref, SNP_in_gene, MUTATION)))
            }
      }
    Out
    }

  # Stop the parallel cluster
  stopCluster(cl)

  Out <- full_join(select(SNP_List, SNP), Temp)

  Out$Gene[is.na(Out$Gene)] <- "Non-gene region"
  Out$SNP_Gene[is.na(Out$SNP_Gene)] <- "Non-gene region"

  Out <- select(Out, SNP, Gene, SNP_Gene)

  Out
}

# # # Example of functionality
# library(TGenGenomicTools)
# library(genbankr)
# library(Biostrings)
# library(tidyverse)
# library(foreach)
# library(parallel)
# library(parallelly)
# library(doParallel)
# # # snp_db <- fasta.to.snps(fasta="/scratch/tporter/RSV_20230402_SequencingMethodCompare_SQL/RSV_A_Fasta/TGen_Clinical_RSVA.fasta",
# # #                         ref_seq="/scratch/tporter/RSV_20230630_Phylogenetic_Analysis/",
# # #                         cores=parallelly::availableCores())
# # # snp_db <- read.csv("/scratch/tporter/RSV_20230908_HomoplasticSNPAnalysis/RSV_A_Homoplastic.csv")
# # # ref_seq="/scratch/tporter/RSV_20230825_ROI_Formatting/RSV_A_NC_038235.1.gb"
# # # snp_db <- dplyr::select(snp_db, SNP)
# # # snps.to.amino(snp_db, ref_seq)
# #
# ref_seq <- "/scratch/mfolkerts/references/H37Rv_NC0009623.gb"
# snp_db <- data_frame("SNP" = c("T1600AG", "C160780G", "C160451A", "T491742C", "T1472517C", "T1472517_", "T1472517AT"))
# cores = 3

