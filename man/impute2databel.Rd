\name{impute2databel}
\alias{impute2databel}
\title{impute2databel}
\usage{impute2databel(genofile, samplefile, outfile, makeprob=TRUE, old=FALSE)}
\description{converts IMPUTE-imputed files to DatABEL (filevector) format}
\details{this function converts IMPUTE-imputed files to DatABEL (filevector) format
containing estimated dosages. 
After conversion, two files (outfile.fvi and outfile.fvd), corresponding 
to single filevector object, will appear on the disk; databel_filtered_R 
object connected to these files will be returned to R.}
\value{databel_filtered_R-class object}
\arguments{\item{genofile}{IMPUTE genotype file name}
\item{samplefile}{IMPUTE sample file name}
\item{outfile}{output file name}
\item{makeprob}{wheather probability-files are also to be arranged}
\item{old}{for developers' use}}
