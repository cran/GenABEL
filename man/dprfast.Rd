\name{dprfast}
\alias{dprfast}
\title{Estimates D' between multiple markers}
\description{
	Given a set of SNPs, computes a matrix of D' 
}
\usage{
dprfast(data, snpsubset, idsubset)
}
\arguments{
  \item{data}{object of \code{\link{snp.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
}
\details{
	The function is based on slightly modified code of Hao et al.
}
\value{
	A (Nsnps X Nsnps) matrix giving D' values
	between a pairs of SNPs above the diagonal and number of SNP genotype  
	measured for both SNPs below the diagonal 
}
\references{
Hao K, Di X, Cawley S. (2006) LdCompare: rapid computation of single- and multiple-marker D' and genetic coverage. Bioinformatics, 23: 252-254.
}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{rhofast}}
}
\examples{
require(GenABEL.data)
data(ge03d2)
# D's using D'fast
a <- dprfast(ge03d2,snps=c(1:10))
\dontrun{
library(genetics)
# D's using package genetics
b <- LD(as.genotype(ge03d2[,1:10]))$"D'"
# see that the D's are not exactly the same
cor(a[upper.tri(a)],b[upper.tri(b)])
plot(a[upper.tri(a)],b[upper.tri(b)])
}
}
\keyword{htest}% at least one, from doc/KEYWORDS
