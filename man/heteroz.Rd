\name{heteroz}
\alias{heteroz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to compute average cross-marker heterozygosity}
\description{
This function computes average heterozygosity within (a set of) 
people, across multiple markers. Can be used for Quality Control
(contamination checks)
}
\usage{
heteroz(x, snpsubset, idsubset, bcast=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of \link{gwaa.data-class} or \link{snp.data-class}}
  \item{snpsubset}{Subset of SNPs to compute average heterozygosity}
  \item{idsubset}{People for whom average heterozygosity should be computed}
  \item{bcast}{Show progress for every \code{bcast} people processed}
}
\details{
}
\value{
  A list with following items
  \item{mean}{Mean heterosygosity}
  \item{sd}{Standard deviation from the mean}
  \item{sem}{Standard error of the mean}
  \item{nmea}{Number of measured SNPs}
  \item{call}{Call rate}
  \item{idsubset}{As in the call}
  \item{snpsubset}{As in the call}
}
\references{}
\author{Yurii Aulchenko}
\note{
}
\seealso{
\code{\link{gwaa.data-class}},
\code{\link{snp.data-class}}
}
\examples{
data(srdta)
h <- heteroz(srdta[,c(1:100)])
plot(h$mean,h$sem)
}
\keyword{htest}% at least one, from doc/KEYWORDS
