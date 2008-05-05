\name{sset}
\alias{sset}
\title{Internal use function for class snp.mx-class}
\description{
Interface to C function sset subsetting genotypes from 
\code{\link{snp.mx-class}}
}
\usage{
sset(data, nsnps, nids, list)
}
\arguments{
  \item{data}{genotypic data in internal format}
  \item{nsnps}{no. snps}
  \item{nids}{no. people}
  \item{list}{something internal...}
}
\details{
	Rather simple function which I wrote before 
	discovering \R's setdiff, etc. functions.
}
\value{
	Sub-set from snp.mx-class object
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{snp.mx-class}} 
}
%\examples{
%}
\keyword{misc}% at least one, from doc/KEYWORDS
