\name{convert.snp.mach}
\alias{convert.snp.mach}
\title{function to convert genotypic data from MACH format to internal data format}
\description{
Converts genotypic data from MACH format to internal genotypic data formated file
}
\usage{
convert.snp.mach(pedfile, mapfile, infofile, outfile, quality = 0.9, ...)
}
\arguments{
  \item{pedfile}{
	File with genotypic data from MACH (geno or mlgeno)
	}
  \item{mapfile}{
	Name of the map file
	}
  \item{infofile}{
	Name MACH info-file
	}
  \item{outfile}{
	Output data file
	}
  \item{quality}{
	Drop the SNPs with quality (as specified in 6th column of info-file) lower than this threshold.
	}
  \item{...}{
	Other arguments passed to \code{\link{convert.snp.ped}}
	}
}
\details{
	This is a simple script converting the MACH data with  \code{\link{convert.snp.ped}},
	re-loading data, and filtering the snp.data object based on quality as specified in 
	MACH info-file
}
\value{
	Does not return any value, but writes file with GenABEL raw data
}
%\references{}
\author{Yurii Aulchenko}
\note{
	The function does not check if "outfile" already exists, thus
	it is always over-written
}
\seealso{
	\code{\link{load.gwaa.data}},
	\code{\link{convert.snp.illumina}},
	\code{\link{convert.snp.text}},
	\code{\link{convert.snp.ped}},
	\code{\link{convert.snp.tped}}
}
\examples{
#
# convert.snp.mach(ped="pedin.18",map="map.18",out="genos.raw")
#
}
\keyword{IO}
