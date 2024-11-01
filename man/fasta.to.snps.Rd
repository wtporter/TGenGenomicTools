% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fasta.to.snps.R
\name{fasta.to.snps}
\alias{fasta.to.snps}
\title{Fasta to SNPS}
\usage{
fasta.to.snps(fasta, ref_seq, cores = parallelly::availableCores())
}
\arguments{
\item{fasta}{A path to the fasta file}

\item{ref_seq}{A path to a reference sequence (.gb) format.}

\item{cores}{Option to}
}
\value{
a data frame with all SNPs within the sequence and sequence name
}
\description{
This function takes a fasta and compares it to reference sequence and computes the SNPs within the sequence. This was built for viral genomes so may struggle with other organisms.
}
