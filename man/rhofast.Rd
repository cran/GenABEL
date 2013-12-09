\name{rhofast}
\alias{rhofast}
\title{Estimates rho between multiple markers}
\description{
	Given a set of SNPs, computes a matrix of rho 
}
\usage{
rhofast(data, snpsubset, idsubset)
}
\arguments{
  \item{data}{object of \code{\link{snp.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
}
\details{
	Rho is the measure of association described by N. Morton and A. Collins
 	(see reference). 
	The function is based on slightly modified code of Hao et al.
}
\value{
	A (Nsnps X Nsnps) matrix giving rho values
	between a pairs of SNPs above the diagonal and 
	Kij below the diagonal.
}
\references{
Collins A, Morton NE. (1998) Mapping a disease locus by allelic association. PNAS, 17:1741-1745.

Hao K, Di X, Cawley S. (2006) LdCompare: rapid computation of single- and multiple-marker rho and genetic coverage. Bioinformatics, 23: 252-254.
}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{r2fast}}
}
\examples{
require(GenABEL.data)
data(ge03d2)
# rhos using rhofast
a <- rhofast(ge03d2,snps=c(1:10))
\dontrun{
library(genetics)
# rhos using package genetics
b <- LD(as.genotype(ge03d2[,1:10]))$"R^2"
# see that the rhos are not exactly the same
cor(a[upper.tri(a)],b[upper.tri(b)])
plot(a[upper.tri(a)],b[upper.tri(b)])
}
}
\keyword{htest}% at least one, from doc/KEYWORDS
