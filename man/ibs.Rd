\name{ibs}
\alias{ibs}
\title{Computes (average) Idenity-by-State for a set of people and markers}
\description{
	Given a set of SNPs, computes a matrix of average IBS for a group of 
	people 
}
\usage{
ibs(data, snpsubset, idsubset, weight="no")
}
\arguments{
  \item{data}{object of \code{\link{snp.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
  \item{weight}{"no" for direct IBS computations, "freq" to weight by allelic frequency}
}
\details{
	This function facilitates quality control of genomic data. 
	E.g. people with exteremly high (close to 1) IBS may indicate 
	duplicated samples (or twins), simply high values of IBS may 
	indicate relatives. 

	When weight "freq" is used, IBS is weighted by (1-freq) where 
	freq is the expected (under HWE) frequency of the genotype,
	if the genotypes are identical, and allele, if alleles are 
	identical.

	ibs() operation may be very lengthy for a large number of people.
}
\value{
	A (Npeople X Npeople) matrix giving average IBS values
	between a pair below the diagonal and number of SNP genotype  
	measured for both members of the pair above the diagonal 
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{check.marker}},
\code{\link{summary.snp.data}},
\code{\link{snp.data-class}}
}
\examples{
data(ge03d2c)
a <- ibs(data=ge03d2c,ids=c(1:10),snps=c(1:1000))
a
# compute IBS based on a random sample of 1000 autosomal marker
a <- ibs(ge03d2c,snps=sample(ge03d2c@gtdata@snpnames[ge03d2c@gtdata@chromosome!="X"],1000,replace=FALSE))
mds <- cmdscale(as.dist(1-a))
plot(mds)
# identify smaller cluster of outliers
km <- kmeans(mds,centers=2,nstart=1000)
cl1 <- names(which(km$cluster==1))
cl2 <- names(which(km$cluster==2))
if (length(cl1) > length(cl2)) cl1 <- cl2;
cl1
# PAINT THE OUTLIERS IN RED
points(mds[cl1,],pch=19,col="red")
}
\keyword{htest}% at least one, from doc/KEYWORDS
