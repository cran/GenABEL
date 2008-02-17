\name{summary.snp.data}
\alias{summary.snp.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to summary GWAA data}
\description{
Provides summary of an object of class \code{\link{snp.data-class}}.
Number of observed genotypes, allelic frequency, genotypic distrbution,
P-value of the exact test for HWE and chromosome are listed
}
\usage{
summary.snp.data(object, ...)
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
	data, but other statistics (freqencies, calls, ...) are based on all data.
	Statistics for Y-chromosome are based on male-only. P-HWE is not defined 
	for mt- and Y- markers (set to 1.0).
}
\seealso{
\code{\link{summary.gwaa.data}},
\code{\link{snp.data-class}}
}
\examples{
data(srdta)
summary(srdta@gtdata[,1:20])
}
\keyword{htest}% at least one, from doc/KEYWORDS

