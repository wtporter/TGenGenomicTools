% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.coding.regions.R
\name{extract.coding.regions}
\alias{extract.coding.regions}
\title{Extract Coding Regions}
\usage{
extract.coding.regions(
  DNAStringSet,
  ref_seq,
  cores = parallelly::availableCores()
)
}
\arguments{
\item{DNAStringSet}{A DNAStringSet object with sequences to be analyzed.}

\item{ref_seq}{A .gb file opened with genbankr::readGenBank to be used as the reference sequence. Check to make sure the @sequence and @transcript slots are full.}

\item{cores}{Option to specify number of cores to use.}
}
\value{
A dataframe with names of the sequences and all coding regions.
}
\description{
This function takes in a DNAStringSet of sequences and uses a reference file to extract the coding region.
}
