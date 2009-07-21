\name{convert.snp.text}
\alias{convert.snp.text}
\title{function to convert integer genotypic data file to raw internal data formated file}
\description{
Converts integer genotypic data file to raw internal data formated file
}
\usage{
convert.snp.text(infile, outfile, bcast = 10000)
}
\arguments{
  \item{infile}{
	Input data file name
	}
  \item{outfile}{
	Output data file
	}
  \item{bcast}{
	Reports progress after reading bcast portion of SNPs
	}
}
\details{
Input genotypic data file contains all kind of genetic information. The first line
of this file contains IDs of all study subjects. The second line gives names of all SNPs 
in the study. The third line list the chromosomes the SNPs belong to. Sequential numbers
are used for autosomes and "X" (capital!) is used for the sex-chromosome. The forth 
line lists genomic position of the SNPs, in order which is the same as order in the line 2. 
The genomic position can be chromosome-specific (each chromosome starts with "0") or, better, 
a true genomic position (chromosome 1 starts with 0 and chromosome 2 continues at the point
chromosome 1 ends). 

Starting with the line five, genetic data are presented. The 5th line contains the data for SNP,
which is listed first on the second line. The first column of this line specifies the 
genotype for the person, who is listed first on the line 1; the second column gives the 
genotype for the second person, so on. The genotypes are coded as 0 (missing), 1 (for AA), 
2 (for AB) and 3 (for BB). Here is a small example:

289982 325286 357273 872422 1005389

SNP-1886933 SNP-2264565 SNP-2305014

1 1 1

825852 2137143 2585920

3 3 3 3 2

3 2 3 3 3

2 2 1 1 1

In this example, we can see that SNP-2305014 (number 3 in the second line) is located on 
chromosome 1 at the position 2585920. If we would like to know what is genotype of 
person with ID 325286 (second in the first line), we need to take second column and 
the third line of the genotypic data. This cell contains 1, thus, person 325286 has 
genotype "AA" at SNP-2305014.

In the event that you do not want to use a map for some reason (such as prior ordering 
of the polymorphisms in the genotype file), make a dummy map-line, which 
contains order information.

The above described genotypic data file is (more or less) human-readable; actually, 
to achieve the aim of effective data storage GWAA package uses internal format. 
In this format, four genotypes are stored in single byte; "raw" data format of R 
is used.
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
	\code{\link{convert.snp.illumina}},
	\code{\link{convert.snp.ped}}, 
	\code{\link{convert.snp.mach}},
	\code{\link{convert.snp.tped}}
}
\examples{
#
# convert.snp.text("genos.dat","genos.raw")
#
}
\keyword{IO}
