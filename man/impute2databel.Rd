\name{impute2databel}
\alias{impute2databel}
\title{converts IMPUTE-imputed files to DatABEL (filevector) format}
\usage{
  impute2databel(genofile, samplefile, outfile,
    makeprob = TRUE, old = FALSE, dataOutType = "FLOAT")
}
\arguments{
  \item{genofile}{IMPUTE genotype file name}

  \item{samplefile}{IMPUTE sample file name}

  \item{outfile}{output file name}

  \item{makeprob}{whether probability-files are also to be
  arranged}

  \item{old}{for developers' use only}

  \item{dataOutType}{the output data type, either "FLOAT"
  or "DOUBLE" (or another DatABEL/filevector type)}
}
\value{
  'databel-class' object
}
\description{
  this function converts IMPUTE-imputed files to DatABEL
  (filevector) format containing estimated dosages. After
  conversion, two files (outfile.fvi and outfile.fvd),
  corresponding to a single filevector object, will appear
  on the disk; a 'databel-class' object connected to these
  files will be returned to R.
}

