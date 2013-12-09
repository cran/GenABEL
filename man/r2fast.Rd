\name{r2fast}
\alias{r2fast}
\title{Estimates r2 between multiple markers}
\description{
	Given a set of SNPs, computes a matrix of r2 
}
\usage{
r2fast(data, snpsubset, cross.snpsubset, idsubset)
}
\arguments{
  \item{data}{object of \code{\link{snp.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{cross.snpsubset}{Parameter allowing parallel implementation. Not to be used normally.
		If supplied together with snpsubset, the r2 for all pairs between 
		snpsubset and cross.snpsubset computed. 
		}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
}
\details{
	The function is based on slightly modified code of Hao et al.
}
\value{
	A (Nsnps X Nsnps) matrix giving r2 values
	between a pairs of SNPs above the diagonal and number of SNP genotype  
	measured for both SNPs below the diagonal 
}
\references{
Hao K, Di X, Cawley S. (2006) LdCompare: rapid computation of single- and multiple-marker r2 and genetic coverage. Bioinformatics, 23: 252-254.
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
# r2s using r2fast
a <- r2fast(ge03d2,snps=c(1:10))
\dontrun{
library(genetics)
# r2s using package genetics
b <- LD(as.genotype(ge03d2[,1:10]))$"R^2"
# see that the r2s are not exactly the same
cor(a[upper.tri(a)],b[upper.tri(b)])
plot(a[upper.tri(a)],b[upper.tri(b)])
}
}
\keyword{htest}% at least one, from doc/KEYWORDS
