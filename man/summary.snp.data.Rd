\name{summary.snp.data}
\alias{summary.snp.data}
\alias{summary.snp.data_old}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to summary GWAA data}
\description{
Provides summary of an object of class \code{\link{snp.data-class}}.
Chromosome, map position, allele coding, number of observed genotypes, 
allelic frequency, genotypic distrbution,
P-value of the exact test for HWE, Fmax (estimate of deviation from HWE, 
allowing meta-analysis) and LRT P-value for HWE test are listed
}
\usage{
	\method{summary}{snp.data}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{snp.data object}
  \item{\dots}{additional arguments (not used)}
}
%\details{}
\value{
	Data frame summary for snp.data object
}
\references{
Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 887-93.
}
\author{Yurii Aulchenko}
\note{
	The P-values reported for X-chromosome are based on analysis of female 
	data, but other statistics (frequencies, calls, ...) are based on all data.
	Statistics for Y-chromosome are based on male-only. P-HWE is not defined 
	for mt- and Y- markers (set to 1.0).
}
\seealso{
\code{\link{summary.gwaa.data}},
\code{\link{snp.data-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
summary(gtdata(srdta[,1:20]))
}
\keyword{htest}% at least one, from doc/KEYWORDS

