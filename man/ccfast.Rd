\name{ccfast}
\alias{ccfast}
\title{fast case-control analysis}
\description{
Fast case-control analysis by computing chi-square test 
from 2x2 (allelic) or 2x3 (genotypic) tables
}
\usage{
ccfast(y, data, snpsubset, idsubset, times=1, quiet=FALSE,bcast=10,
		clambda=TRUE,propPs=1.0)
}
\arguments{
  \item{y}{character name of the vector of case-control status. Cases are denoted as 1 and controls as 0.}
  \item{data}{An object of \code{\link{gwaa.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
  \item{times}{If more then one, the number of replicas to be used in derivation of 
		empirical genome-wide significance. See \code{\link{emp.qtscore}}, which
		calls qtscore with times>1 for details}
  \item{quiet}{do not print warning messages}
  \item{bcast}{If the argument times > 1, progress is reported once in bcast replicas}
  \item{clambda}{If inflation facot Lambda is estimated as lower then one, this parameter 
		controls if the original P1df (clambda=TRUE) to be reported in Pc1df, 
		or the original 1df statistics is to be multiplied onto this "deflation" 
		factor (clambda=FALSE). 
		If a numeric value is provided, it is used as a correction factor.}
  \item{propPs}{proportion of non-corrected P-values used to estimate the inflation factor Lambda,
		passed directly to the \code{\link{estlambda}}}
}
%\details{
%}
\value{
  Object of class \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{emp.ccfast}},
\code{\link{plot.scan.gwaa}},
\code{\link{scan.gwaa-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
a <- ccfast("bt",data=srdta,snps=c(1:10),ids=c(1:100))
a
a <- ccfast("bt",data=srdta)
plot(a)
}
\keyword{htest}% at least one, from doc/KEYWORDS
