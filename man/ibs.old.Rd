\name{ibs.old}
\alias{ibs.old}
\title{Computes (average) Idenity-by-State for a set of people and markers}
\description{
	Given a set of SNPs, computes a matrix of average IBS for a group of 
	people 
}
\usage{
ibs.old(data, snpsubset, idsubset, weight="no")
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

	When weight "freq" is used, IBS for a pair of people i and j is computed as

	\deqn{
	f_{i,j} = \Sigma_k \frac{(x_{i,k} - p_k) * (x_{j,k} - p_k)}{(p_k * (1 - p_k))}
	}

	where k changes from 1 to N = number of SNPs GW, \eqn{x_{i,k}} is 
	a genotype of ith person at the kth SNP, coded as 0, 1/2, 1 and 
	\eqn{p_k} is the frequency 
	of the "+" allele. This apparently provides an unbiased estimate of the 
	kinship coefficient.

	Only with "freq" option monomorphic SNPs are regarded as non-informative.

	ibs() operation may be very lengthy for a large number of people.
}
\value{
	A (Npeople X Npeople) matrix giving average IBS (kinship) values
	between a pair below the diagonal and number of SNP genotype  
	measured for both members of the pair above the diagonal. 

	On the diagonal, homozygosity (0.5+inbreeding) is provided.
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
require(GenABEL.data)
data(ge03d2c)
a <- ibs(data=ge03d2c,ids=c(1:10),snps=c(1:1000))
a
# compute IBS based on a random sample of 1000 autosomal marker
a <- ibs(ge03d2c,snps=sample(ge03d2c@gtdata@snpnames[ge03d2c@gtdata@chromosome!="X"],
			1000,replace=FALSE),weight="freq")
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
