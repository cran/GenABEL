\name{impute2mach}
\alias{impute2mach}
\title{converts IMPUTE to MACH files}
\usage{
  impute2mach(genofile, infofile, samplefile, machbasename,
    maketextdosefile = TRUE, ...)
}
\arguments{
  \item{genofile}{IMPUTE genotype file name}

  \item{infofile}{IMPUTE info file name}

  \item{samplefile}{IMPUTE sample file name}

  \item{machbasename}{base name for MACH-formatted outputs}

  \item{maketextdosefile}{whether text dosefile is to be
  generated (if not, only filevector (*.fvi / *.fvd) files,
  usable with ProbABEL, will be generated)}

  \item{...}{arguments passed to
  \link{extract.annotation.impute} (DO CHECK documentation,
  otherwise your annotation may be skrewed up!)}
}
\value{
  nothing returned except files generated on the disk
}
\description{
  function to convert IMPUTE files to MACH format
}
\author{
  Yurii Aulchenko
}
\keyword{IO}
\keyword{manip}

