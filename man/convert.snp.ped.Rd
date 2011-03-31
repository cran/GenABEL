\name{convert.snp.ped}
\alias{convert.snp.ped}
\title{function to convert genotypic data in pedigree fromat (+map) to internal data format}
\description{
Converts genotypic data in a variety of pedigree fromats (+map) to internal genotypic data formated file
}
\usage{
convert.snp.ped(pedfile, mapfile, outfile, format = "premakeped", traits = 1, 
	strand = "u", bcast = 10000000, wslash=FALSE, mapHasHeaderLine=TRUE)
}
\arguments{
  \item{pedfile}{
	Pre-makeped linkage genotypic data file name (no header line)
	}
  \item{mapfile}{
	Name of the map file (note that by default header line should be included, 
	see 'mapHasHeaderLine' to modify that)
	}
  \item{outfile}{
	Output data file
	}
  \item{format}{
	Input data format, either "premakeped" (default, also works with Merlin 
	files), or "mach"
	}
  \item{traits}{
	How many traits are specified in the pedigree file (usually 1 -- affection
	-- or 2 -- affection and liability).
	Has no effect when format = "mach".
	}
  \item{strand}{
	Specification of strand, one of "u" (unknown), "+", "-" or "file".
	In the latter case, map-file should contain an extended map (the one 
	including strand and coding)
	}
  \item{bcast}{
	Reports progress after reading bcast portion of SNPs
	}
  \item{wslash}{
	Whether alleles are separated with slash (is true for Mach/Merlin format), 
	otherwise it is assumed that alleles are separated with space. When 
	wslash=T it is assumed that genotypes are coded with single characters, 
	separated with slash (no spaces), e.g. "A/G", and not "A/ G" or "A / G".
	}
  \item{mapHasHeaderLine}{
	Whether map-file has a header line
	}
}
\details{
	Pedfile must be standard pre-makeped/Merlin linkage file, or a Mach file. 
	In pre-makeped linkage file, columns are
	
	ped id fa mo sex trait snp1.allele1 snp1.allele2 snp2.allele1 snp2.allele2 ...

	For example

	1 1 0 0 1 2 A A G T ...

	1 2 0 0 1 0 A G T T ...

	1 3 0 0 2 1 A A T T ...

	...

	Would imply that persons 1, 2 and 3 are "founders" (which would be typical for a 
	case-control study), 1 and 2 are males and 3 is female. Person 1 is homozygous for 
	allele 1 at locus 1 and heterozygous at locus 2. Person 2 is heterozygous at both loci. 
	Person 3 is homozygous for allele 2 at locus 1 and allele 1 at locus 2.

	Only the second and the marker columns are used, thus make sure the IDs are unique!

	Accepted allele codes: 1/2, A/B, A/T, A/G, A/C, T/G, T/C, G/C, A/-, T/-, G/-, C/-.
	Here, "-" stands of a deletion. 

	The map file is standard Merlin map. For example:

	chrom	name	position

	18	rs679153	2859916	

	18	rs9965482	2860891

	Says that locus 1 is named rs679153 and located at chromosome 18 position 2859916.
	Locus 2 (rs9965482) is located at chromosome 18, position 2860891.

	In extended map format, there should be 4th column specifying the strand

	chrom name       position strand

	18    rs679153   2859916   -

	18    rs9965482  2860891   +

	Accepted strand coding: +, -, u (unknown)

	Please note that by default the header line (e.g. "chrom	name	position") SHOULD 
	be present in your file, though you can use 'mapHasHeaderLine' argument to modify 
	this behavior
	
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
	\code{\link{convert.snp.mach}},
	\code{\link{convert.snp.text}},
	\code{\link{convert.snp.tped}}
}
\examples{
#
# convert.snp.ped(ped="pedin.18",mapfile="map.18",out="genos.raw")
#
}
\keyword{IO}
