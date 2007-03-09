\name{convert.snp.ped}
\alias{convert.snp.ped}
\title{function to convert genotypic data in pre-makeped linkage fromat (+map) to internal genotypic data formated file}
\description{
Converts genotypic data in pre-makeped linkage fromat (+map) to internal genotypic data formated file
}
\usage{
convert.snp.ped(pedfile, mapfile, outfile, bcast = 10000)
}
\arguments{
  \item{pedfile}{
	Pre-makeped linkage genotypic data file name
	}
  \item{mapfile}{
	Mega2 map file
	}
  \item{outfile}{
	Output data file
	}
  \item{bcast}{
	Reports progress after reading bcast portion of SNPs
	}
}
\details{
	Pedfile must be standard pre-makeped linkage file. In this file,
	columns are
	
	ped id fa mo sex affection

	Sex is coded as 1=male and 2=female. Affection status is not used. For example

	1 1 0 0 1 2 1 1 1 2 

	1 2 0 0 1 0 1 2 1 2

	1 3 0 0 2 1 2 2 1 1

	Would imply that persons 1, 2 and 3 are "founders" (which would be typical for a 
	case-control study), 1 and 2 are males and 3 is female. Person 1 is homozygous for 
	allele 1 at locus 1 and heterozygous at locus 2. Person 2 is heterozygous at both loci. 
	Person 3 is homozygous for allele 2 at locus 1 and allele 1 at locus 2.

	The map file is standard Mega2 map. For example:

	chrom	kosambi	name

	18	2859916	rs679153

	18	2860891	rs9965482

	Says that locus 1 is named rs679153 and located at chromosome 18 position 2859916.
	Locus 2 (rs9965482) is located at chromosome 18, position 2860891.
	
}
\value{
	Does not return any value
}
%\references{}
\author{Yurii Aulchenko}
\note{
	The function does not check if "outfile" already exists, thus
	it is always over-written
}
\seealso{
	\code{\link{load.gwaa.data}},
	\code{\link{convert.snp.text}}
}
\examples{
#
# convert.snp.ped(ped="pedin.18",map="map.18",out="genos.raw")
#
}
\keyword{IO}
