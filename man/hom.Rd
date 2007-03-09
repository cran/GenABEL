\name{hom}
\alias{hom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to compute average homozygosity within a person}
\description{
This function computes average homozygosity for a set of 
people, across multiple markers. Can be used for Quality Control
(e.g. contamination checks)
}
\usage{
	hom(data, snpsubset, idsubset, weight="no")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Object of \link{gwaa.data-class} or \link{snp.data-class}}
  \item{snpsubset}{Subset of SNPs to be used}
  \item{idsubset}{People for whom average homozygosity is to be computed}
  \item{weight}{When "no", homozygosity is computed as a proportion of 
		homozygous genotypes. When "freq", homozygosity is computed 
		as a sum over one minus the frequency of the homozygous allele, 
		and then divided by the number of genotypes.}
}
\details{
}
\value{
	A matrix with rows corresponding to the ID names and colums
	showing the number of genotypes measured (NoMeasured) and 
	homozygosity (Hom).
}
\references{}
\author{Yurii Aulchenko}
\note{
}
\seealso{
\code{\link{ibs}},
\code{\link{gwaa.data-class}},
\code{\link{snp.data-class}}
}
\examples{
data(ge03d2)
h <- hom(ge03d2[,c(1:100)])
homsem <- h[,"Hom"]*(1-h[,"Hom"])/h[,"NoMeasured"]
plot(h[,"Hom"],homsem)
}
\keyword{htest}% at least one, from doc/KEYWORDS
