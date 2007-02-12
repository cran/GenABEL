\name{check.marker}
\alias{check.marker}
\title{function to do genotypic quality control}
\description{
This function helps selecting the marker which should enter into 
GWA analysis based on call rate, minor allele frequency, 
value of the chi-square test for Hardy-Weinberg equilibrium, 
and redudndance, defined as concordance between the distributions 
of the genotypes (including missing values).
}
\usage{
check.marker(data, snpsubset, idsubset, callrate = 0.95, 
		maf = 0.01, p.level = -1, fdrate = 0.2, 
		hweidsubset, redundant = "bychrom", 
		minconcordance = 2, qoption = "bh95", 
		bcast = 10000)
}
\arguments{
  \item{data}{gwaa.data or snp.data object}
  \item{snpsubset}{a subset of SNPs to check (names, indexes, logical), default is all from \code{data}}
  \item{idsubset}{a subset of people to check (names, indexes, logical), default is all from \code{data}}
  \item{callrate}{cut-off call rate}
  \item{maf}{cut-off Minor Allele Frequency}
  \item{p.level}{cut-off p-value in check for Hardy-Weinberg Equilibrium. If negative, FDR is applied}
  \item{fdrate}{cut-off FDR level in check for Hardy-Weinberg Equilibrium}
  \item{hweidsubset}{a subset of people to check (names, indexes, logical) to use for HWE check}
  \item{redundant}{if "bychrom", redundancy is checked within chromosomes; "all" -- all pairs of markers; "no" -- no redundancy checks}
  \item{minconcordance}{a parameter passed to "redundant" function. If "minconcordance" is > 1.0 
			only pairs of markers which are exactly the same, including NA pattern,
			are considered as redundant; if 0 < "minconcordance" < 1, then pairs 
			of markers with concordance > "minconcordance" are considered redundant.
			see \code{\link{redundant}} for details. Note that if "minconcordance" < 1
			the program will take much longer time to run}
  \item{qoption}{if "bh95", BH95 FDR used; if "storey", qvalue package is used}
  \item{bcast}{report progress every bcast SNPs}
}
%\details{
%}
\value{
  Object of class \code{\link{check.marker-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
	\code{\link{check.trait}},
	\code{\link{plot.check.marker}},
	\code{\link{summary.check.marker}},
	\code{\link{redundant}},
	\code{\link{HWE.show}},
	\code{\link{check.marker-class}}
}
\examples{
data(srdta)
mc <- check.marker(data=srdta)
names(mc)
mc$nohwe
mc <- check.marker(data=srdta@gtdata[,1:100],call=0.98,maf=0.01,minconcordance=0.9,fdr=.1)
summary(mc)
HWE.show(data=srdta,snps=mc$nohwe)
plot(mc)
}
\keyword{misc}% at least one, from doc/KEYWORDS
