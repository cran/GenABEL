\name{convert.snp.tped}
\alias{convert.snp.tped}
\title{function to convert genotypic data in transposed-ped format
  (.tped and .tfam) to internal genotypic data formatted file}
\description{
Converts genotypic data in transposed-ped format (.tped and .tfam) to internal genotypic data formatted file
}
\usage{
convert.snp.tped(tpedfile, tfamfile, outfile,strand = "u", bcast = 10000)
}
\arguments{
  \item{tpedfile}{
	Name of transposed-ped format (.tped) file to read
	}
  \item{tfamfile}{
	Name of individual data (.tfam) file to read
	}
  \item{outfile}{
	Name for output data file
	}
  \item{strand}{
	Specification of strand, one of "u" (unknown), "+", "-" or "file".
	In the latter case, extra column specifying the strand (again, one of 
	"u", "+", or "-") should be included on the tpedfile.
	}
  \item{bcast}{
	Reports progress every time this number of SNPs have been read
	}
}
\details{
  The transposed-ped file format may be preferred when extremely large numbers of
  markers have been genotyped.  This file format is supported by plink!
  See http://pngu.mgh.harvard.edu/~purcell/plink/ for details.

  The conversion is performed by C++ code that is both fast and memory
  efficient.
  
  The genotype data are stored in the main transposed-ped format file,
  usually with a .tped file extension.  If there are NSNP markers
  genotyped in NIND individuals, this file has NSNP rows and 4+NIND*2
  columns.  There is one row per marker, and no header.  The first four columns are:

  Chromosome

  Marker name (e.g. rs number)

  Genetic position (in Morgans)

  Physical position (in bp)
  
  These are followed by two columns per individual, which contain the
  genotype, coded as two characters.  The `0' character is used for
  missing data.  For example, a file containing data for six individuals
  genotyped at two SNPs would look like:
  
  1 rs1234 0 5000650 A A 0 0 C C A C C C C C

  1 rs5678 0 5000830 G T G T G G T T G T T T  

  In this example, the second individual is missing data for SNP rs1234,
  etc.  The alleles can be coded by any two distinct characters,
  e.g. 'C' and 'G', or '1' and '2'.  The '0' character is reserved for missing data,
  and each individual genotype must be either complete, or completely
  missing.  In the current implementation, only the physical positions
  of the SNPs are read, and the genetic positions are ignored.
  
  The indices for the columns are stored in a separate file, usually
  with a .tfam file extension.  Traditionally, this file has six
  columns, and no header.  In the current implementation, only the
  second column is used.  This column must contain the individual id.  Other
  columns are ignored.
	
}
\value{
	Does not return any value
}
%\references{}
\author{Toby Johnson <toby.johnson@unil.ch>}
\note{
	The function does not check if "outfile" already exists, thus
	it is always over-written
}
\seealso{
	\code{\link{convert.snp.ped}}, 
	\code{\link{convert.snp.illumina}},
	\code{\link{convert.snp.text}},
	\code{\link{convert.snp.mach}},
	\code{\link{load.gwaa.data}}
}
\examples{
#
# convert.snp.tped("c21.tped",map="c21.tfam",out="c21.raw")
#
}
\keyword{IO}
