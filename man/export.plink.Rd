\name{export.plink}
\alias{export.plink}
\title{Export GenABEL data in PLINK format}
\usage{
  export.plink(data, filebasename = "plink",
    phenotypes = "all", transpose = TRUE,
    export012na = FALSE, ...)
}
\arguments{
  \item{data}{GenABEL data object of 'gwaa.data'-class to
  be exported.}

  \item{filebasename}{base file name for exported data,
  extensions '.ped', '.map' and '.phe' (for phenotype file)
  are added for specific output files.}

  \item{phenotypes}{NULL (no phenotypes exported), "all"
  (default) for all phenotypes or a vector of character
  with names of phenotypes to be exported.}

  \item{transpose}{if TRUE (default), 'tped' files will be
  produced, else 'ped' files are produced.}

  \item{export012na}{if TRUE, export in numeric (0, 1, 2,
  NA) format, as opposed to ATGC format (default: FALSE).}

  \item{...}{arguments passed to
  \code{\link{export.merlin}}.}
}
\description{
  Export GenABEL data in PLINK format. This function is a
  simple wrapper to the \code{\link{export.merlin}}
  function with specific arguments + few lines of code to
  export phenotypes
}
\author{
  Yurii Aulchenko
}
\keyword{IO}

