\name{mach2databel}
\alias{mach2databel}
\title{mach2databel}
\usage{mach2databel(imputedgenofile, mlinfofile, outfile, isprobfile=FALSE)}
\description{converts MACH-imputed files to DatABEL (filevector) format}
\details{This function converts mach-imputed files to \code{DatABEL} (filevector) format.
After conversion, two files (outfile.fvi and outfile.fvd), corresponding 
to single filevector object, will appear on the disk; 'databel-class' 
object connected to these files will be returned to R}
\value{databel-class object}
\author{Yurii Aulchenko}
\keyword{IO}
\keyword{manip}
\arguments{\item{imputedgenofile}{MACH mldose (or mlprob) file name}
\item{mlinfofile}{MACH mlinfo file name}
\item{outfile}{output file name}
\item{isprobfile}{whether imputedgenofile is a prob-file 
(default is FALSE, that is dose-file assumed)}}
