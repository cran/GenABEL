\name{summary.gwaa.data}
\alias{summary.gwaa.data}
\title{function to summarise GWAA data}
\description{
Summary of phenotypic and genotypic parts of GWAA data
}
\usage{
\method{summary}{gwaa.data}(object, ...) 
}
\arguments{
  \item{object}{object of class \code{\link{gwaa.data-class}}}
  \item{\dots}{additional arguments (not used)}
}
%\details{}
\value{
	Returns list with two elements: 
	\item{pheno}{Summary for phenotypic part of gwaa.data object}
	\item{geno}{Summary for genotypic part of gwaa.data object}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{\code{\link{summary.snp.data}}}
\examples{
require(GenABEL.data)
data(srdta)
# be prepared : long output!
summary(srdta)
}
\keyword{htest}% at least one, from doc/KEYWORDS

