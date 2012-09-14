\name{arrange_probabel_phe}
\alias{arrange_probabel_phe}
\title{arranges ProbABEL phenotype-file}
\usage{
  arrange_probabel_phe(modelterms, phedata, gendata,
    file = "probabel.PHE")
}
\arguments{
  \item{modelterms}{vector of character, which specifies
  the variables to be included into ProbABEL
  phenotype-file. Should contain, and start with 'id'
  column, which should provide the same ID codes as these
  in gendata}

  \item{phedata}{phenotypic data (matrix, data.frame, or
  \link{gwaa.data-class} object)}

  \item{gendata}{genetic data to be used with ProbABEL,
  either databel-class object or name of the (index/data)
  file containing filevector data for ProbABEL}

  \item{file}{name of the ProbABEL phenotype file}
}
\value{
  file with phenotypes ready for use with ProbABEL
}
\description{
  Function to arrange ProbABEL phenotype-file; it takes
  phenotypic data as input and aligns that with genotypic
  data of ProbABEL
}
\author{
  Yurii Aulchenko
}
\keyword{IO}
\keyword{manip}

