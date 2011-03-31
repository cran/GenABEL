\name{convert.snp.illumina}
\alias{convert.snp.illumina}
\title{function to convert genotypic data from Illumina/Affymetrix to internal format}
\description{
Converts genotypic data from Illumina/Affymetrix-like format to internal genotypic data formated file
}
\usage{
convert.snp.illumina(infile, outfile, strand = "+", bcast = 10000000)
}
\arguments{
  \item{infile}{
	Map + genotypic data file name
	}
  \item{outfile}{
	Output data file
	}
  \item{strand}{
	Specification of strand, one of "u" (unknown), "+", "-" or "file".
	In the latter case, extra column specifying the strand (again, one of 
	"u", "+", or "-") should be included on the infile.
	}
  \item{bcast}{
	Reports progress after reading bcast portion of SNP genotypes
	}
}
\details{
	Input file is the one which could be typically obtained from Illumina 
	BeadStudio software. For example:

	Name	Chr	Pos	id1	id2	id3

	rs1001	2	12897	AC	AA	AA

	rs2401	3	12357	AG	GG	AG

	rs123	3	5327	TC	TT	CC

	Here, every row corresponds to a SNP, and each column, starting with the 
	4th, corresponds to a person. 

	When strand information is available (option strand="file"), the file should
	look like

	Accepted allele codes: 1/2, A/B, A/T, A/G, A/C, T/G, T/C, G/C, A/-, T/-, G/-, C/-.
	Here, "-" stands of a deletion. Missing data can be coded as "--" or "00". Make sure 
	that the coding for missing is "00" if you use one of the codings A/-, T/-, G/-, C/-!

	Name	Chr	Pos	Strand	id1	id2	id3

	rs1001	2	12897	+	AC	AA	AA

	rs2401	3	12357	-	AG	GG	AG

	rs123	3	5327	+	TC	TT	CC

	Accepted strand coding: +, -, u (unknown)	

	The procedure always codes genotypes that "0", "1" and "2" correspond to 
	AA, AB, and BB, where B is the less frequent allele. Thus GWA analysis 
	procedures will return effect of the minor allele.
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
	\code{\link{convert.snp.text}},
	\code{\link{convert.snp.mach}},
	\code{\link{convert.snp.tped}}
}
\examples{
#
# convert.snp.illumina(infile="pedin.18",out="genos.raw",strand="+")
#
}
\keyword{IO}
