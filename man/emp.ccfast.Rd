\name{emp.ccfast}
\alias{emp.ccfast}
\title{Genome-wide significance for a case-control GWA scan}
\description{
Genome-wide significance for a case-control GWA scan. Analysis 
function is \code{\link{ccfast}}. 
}
\usage{
emp.ccfast(y, data, snpsubset, idsubset, times = 200, quiet=FALSE,
		bcast = 10)
}
\arguments{

  All arguments are the same as in and passed intact to the \code{\link{ccfast}}.
  See help for this function. 
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
}
\details{
In the analysis of empirical significance, first time the function 
\code{\link{ccfast}} is called and result object is 
saved. Later, the function \code{\link{ccfast}} is called 
\code{times} times with \code{replace=FALSE} in order to generate 
the distribution under the null. Each call, minimal P-value is extracted 
and compared with original P-values. For a particular SNP, empirical 
P-value is obtained as a proportion of times minimal Ps from resampled data 
was less then the original P.

The list elements effB, effAB and effBB are the ones obtained from the 
analysis of the original (not permuted) data set
}
\value{
	Object of class \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{ccfast}},
\code{\link{emp.qtscore}},
\code{\link{scan.gwaa-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
a<-ccfast("bt",data=srdta,snps=c(500:800))
plot(a)
# this does not make sense, as the whole experiment must be analysed, not a small region!
# also, times = 20 is way too small (should be at least 200)
b<-emp.ccfast("bt",data=srdta,snps=c(500:800),bcast=10, times = 20)
plot(b)
# compare qvalues and empirical P
qv<-qvaluebh95(a[,"P1df"])$qval
qv
b[,"P1df"]
plot(qv,b[,"P1df"],xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)
}
\keyword{htest}% at least one, from doc/KEYWORDS
