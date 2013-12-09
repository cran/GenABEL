\name{redundant}
\alias{redundant}
\title{function to do redundancy check}
\description{
Checks marker redundancy, understood as comcordance between 
genotypic distributions (including missing values)
}
\usage{
redundant(data, pairs = "bychrom", minconcordance = 2.0)
}
\arguments{
  \item{data}{gwaa.data or snp.data object}
  \item{pairs}{"bychrom" or "all" to check pairs within chromosome only or genome-wide}
  \item{minconcordance}{find "redundant" pairs of markers with concordance >= "minconcordance".
			If "minconcordance" is more then 1.0, only pairs of markers which 
			are exactly the same (independent of coding), including NA pattern, are considered as 
			redundant. If "minconcordance" is <= 1, the concordance rate is 
			computed as percent of genotypes which are the same, including  
			the genotypes with NA. I.e. if both genotypes are NA, this is 
			counted as a match, if one is NA and other is measured, this is 
			counted as dismatch. Note that option with "minconcordance" <= 1
			takes much longer time to run.}
}
%\details{}
\value{
 A list containing reference SNP as a name and all SNPs which has "the same" 
 genotypic distribution as values: 
  \item{"refSNP1"}{SNP11, SNP12, ...}
  \item{"refSNP2"}{SNP21, SNP22, ...}
  \item{\dots}{etc.}
  \item{"refSNPlast"}{SNPlast1, SNPlast2, ...}
  \item{"all"}{list of all redundant SNPs, which can be dropped from consideration}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{check.marker}}
}
\examples{
require(GenABEL.data)
data(srdta)
redundant(srdta@gtdata)
redundant(srdta@gtdata[,1:50],minconcordance=0.8)
}
\keyword{htest}
