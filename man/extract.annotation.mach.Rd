\name{extract.annotation.mach}
\alias{extract.annotation.mach}
\title{extracts SNP annotation from MACH/HapMap legend files}
\usage{
  extract.annotation.mach(infofile, legendfile,
    chromosome = NA)
}
\arguments{
  \item{infofile}{MACH (ml)info-file name}

  \item{legendfile}{HapMap legend file name}

  \item{chromosome}{chromosome}
}
\value{
  data frame containing annotaton
}
\description{
  This function extracts SNP annotation from MACH info and
  HapMap legend files.
}
\author{
  Yurii Aulchenko
}
\keyword{IO}
\keyword{manip}

