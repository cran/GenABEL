\name{perid.summary}
\alias{perid.summary}
\title{Summary of marker data per person}
\description{
	Produces call rate and heterozygosity per person
}
\usage{
perid.summary(data, snpsubset, idsubset, ... )
}
\arguments{
  \item{data}{object of \code{\link{snp.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
  \item{...}{additional parameters to be passed to \code{\link{hom}}}
}
\details{
	This function facilitates quality control of genomic data. 
	E.g. extreme outliers for heterozygosity indicate possibly 
	contaminated DNA samples, while low call rate of a person 
	may indicate poor DNA quality.
}
\value{
	A matrix, giving per person (row) its' average heterozygosity 
	("Het" column) and call rate ("CallPP"), over all SNPs 
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{check.marker}},
\code{\link{summary.snp.data}},
\code{\link{snp.data-class}}
}
\examples{
require(GenABEL.data)
data(ge03d2c)
a <- perid.summary(data=ge03d2c,snps=c(1:100),ids=c(1:10))
a
a <- perid.summary(data=ge03d2c)
hist(a[,"CallPP"])
hist(a[,"Het"])
}
\keyword{htest}% at least one, from doc/KEYWORDS
