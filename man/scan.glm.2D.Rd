\name{scan.glm.2D}
\alias{scan.glm.2D}
\title{Scans regional data allowing for gene-gene interaction using glm}
\description{
Scans regional data allowing for gene-gene interaction using glm
}
\usage{
scan.glm.2D(formula, family = gaussian(), data, snpsubset, idsubset, 
		bcast = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{character string containing formula to be used in \code{\link[stats]{glm}}. 
		You should put CRSNP argument in the formula,
		to arrange how the SNP from the list would be treated. This allows to 
		put in an interaction term.}
  \item{family}{family to be passed to \code{\link[stats]{glm}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{data}{object of class "gwaa.data"}
  \item{bcast}{show progress every \code{bcast} SNPs}
}
\details{
}
\value{
  Object of class \code{\link{scan.gwaa.2D-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{scan.gwaa.2D-class}},
	\code{\link{scan.haplo.2D}}
}
\examples{
data(srdta)
a <- scan.glm.2D("bt~sex+age+CRSNP",family=binomial(),data=srdta,snps=(1:10),bcast=2)
plot(a)
}
\keyword{htest}% at least one, from doc/KEYWORDS
