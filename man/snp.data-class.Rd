\name{snp.data-class}
\docType{class}
\alias{snp.data-class}
\alias{[,snp.data,ANY,ANY,ANY-method}
\alias{coerce,snp.data,numeric-method}
\alias{coerce,snp.data,character-method}
\alias{coerce,snp.data,genotype-method}
\alias{coerce,snp.data,hsgeno-method}
\alias{show,snp.data-method}
\alias{summary,snp.data-method}

\title{Class "snp.data"}
\description{This class contains objects holding large arrays of single
  nucleotide polymorphism (SNP) genotypes}
\section{Slots}{
	 \describe{
    \item{\code{nbytes}:}{number of bytes used to store data on a SNP}
    \item{\code{nids}:}{number of people}
    \item{\code{male}:}{male code}
    \item{\code{idnames}:}{ID names}
    \item{\code{nsnps}:}{number of SNPs}
    \item{\code{nsnpnames}:}{list of SNP names}
    \item{\code{chromosome}:}{list chromosomes corresponding to SNPs}
    \item{\code{coding}:}{list of nucleotide coding for the SNPs}
    \item{\code{strand}:}{strands of the SNPs}
    \item{\code{map}:}{list SNPs' positions}
    \item{\code{gtps}:}{\code{\link{snp.mx-class}} object used to store genotypes}
	}
}
%\section{Extends}{}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "snp.data", i = "ANY", j = "ANY", drop = "ANY")}: subset
      operations. x[i,j] will select people listed in i and SNPs listed in j.}
    \item{coerce}{\code{signature(from = "snp.data", to = "numeric")}:
      map to codes 0, 1, 2, or NA}
    \item{coerce}{\code{signature(from = "snp.data", to =
	"character")}: map to actual nucleotide codes, e.g. "A/A", "A/G", "G/G", ""}
    \item{coerce}{\code{signature(from = "snp.data", to =
	"genotype")}: map to data frame with \code{\link[genetics]{genotype}}-class
	data, for later use with package \code{genetics}}
    \item{coerce}{\code{signature(from = "snp.data", to =
	"hsgeno")}: map to data frame with allelic data frame, 
	for later use with package \code{haplo.stats}}
    \item{show}{\code{signature(object = "snp.data")}: shows the object.
      Take care that the objects are usually very large!}
    \item{summary}{\code{signature(object = "snp.data")}: calculate
      allele frequencies, genotype frequencies, and chi-square tests for
    Hardy-Weinberg equilibrium. Results are returned as a dataframe}
	 }
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{gwaa.data-class}},
	\code{\link{snp.data}},
	\code{\link{snp.mx-class}}
}
\examples{
data(srdta)
class(srdta)
x <- srdta@gtdata
class(x)
x@nids
x@nsnps
x@idnames[1:12]
x@male[1:12]
x@male[c("p1","p2","p3","p4")]
x@snpnames[1:4]
x@chromosome[1:4]
x@map[1:4]
n4 <- c("rs18","rs655")
n4
x@map[n4]
n4 <- c("rs18","rs65")
n4
x@map[n4]
x@chromosome[n4]
x[1:12,1:4]
summary(x[,1:10])
as.numeric(x[1:12,1:4])
as.numeric(x[c("p1","p3","p4"),c("rs18","rs65")])
as.character(x[c("p1","p3","p4"),c("rs18","rs65")])
as.genotype(x[c("p1","p3","p4"),c("rs18","rs65")])
as.hsgeno(x[c("p1","p3","p4"),c("rs18","rs65")])
}
\keyword{classes}

