\name{mach2databel}
\alias{mach2databel}
\title{converts MACH-imputed files to DatABEL (filevector) format}
\usage{
  mach2databel(imputedgenofile, mlinfofile, outfile,
    isprobfile = FALSE, dataOutType = "FLOAT")
}
\arguments{
  \item{imputedgenofile}{MACH mldose (or mlprob) file name}

  \item{mlinfofile}{MACH mlinfo file name}

  \item{outfile}{output file name}

  \item{isprobfile}{whether imputedgenofile is a prob-file
  (default is FALSE, that is dose-file assumed)}

  \item{dataOutType}{the output data type, either "FLOAT"
  or "DOUBLE" (or other DatABEL/filevector type)}
}
\value{
  databel-class object
}
\description{
  This function converts mach-imputed files to
  \code{DatABEL} (filevector) format. After conversion, two
  files (outfile.fvi and outfile.fvd), corresponding to
  single filevector object, will appear on the disk;
  'databel-class' object connected to these files will be
  returned to R
}
\author{
  Yurii Aulchenko
}
\keyword{IO}
\keyword{manip}

