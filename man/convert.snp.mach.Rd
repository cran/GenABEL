\name{convert.snp.mach}
\alias{convert.snp.mach}
\title{function to convert genotypic data from MACH format to internal data format}
\description{
Converts genotypic data from MACH format to internal genotypic data formated file 
(NOT recommended)
}
\usage{
convert.snp.mach(pedfile, mapfile, infofile, outfile, quality = 0.9, 
					column.quality = 7, strand = "+", ...)
}
\arguments{
  \item{pedfile}{
	File with genotypic data from MACH (geno or mlgeno)
	}
  \item{mapfile}{
	Name of the map file (note that header line should be included)
	}
  \item{infofile}{
	Name MACH info-file
	}
  \item{outfile}{
	Output data file
	}
  \item{quality}{
	Drop the SNPs with quality (as specified in some column of info-file) lower than this threshold.
	}
  \item{column.quality}{
	What column of the info-file provides "quality". Default = 7 os r2; possible values 
	include 6 (average postrior probability).
	}
  \item{strand}{
	Specification of strand, one of "u" (unknown), "+", "-" or "file".
	In the latter case, map-file should contain an extended map (the one 
	including strand and coding). See options to \code{\link{convert.snp.ped}}
	for details.
	}
  \item{...}{
	Other arguments passed to \code{\link{convert.snp.ped}}
	}
}
\details{
	This is a simple script converting the MACH data with  \code{\link{convert.snp.ped}},
	re-loading data, and filtering the snp.data object based on quality as specified in 
	MACH info-file
	
	Note that it is NOT recommended to process imputations results with this procedure, 
	as uncertainty inherent to imputations is lost after processing ("hard" calls are made)!
	For the purposes of GWA analysis, we recommend procedures, which make direct use of 
	estimated allele dose or genotypic probabilities.
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
