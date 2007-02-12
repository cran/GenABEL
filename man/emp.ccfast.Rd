\name{emp.ccfast}
\alias{emp.ccfast}
\title{Bootstrap analysis of case-control GWA scan}
\description{
Bootstrap analysis of case-control GWA scan. Analysis 
function is \code{\link{ccfast}}. 
}
\usage{
emp.ccfast(y, data, snpsubset, idsubset, times = 100, 
		bcast = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  Only last two arguments are specific for \code{emp.ccfast} and others are 
  passed intact to \code{\link{ccfast}}.

  \item{y}{vector of case-control status. Cases are denoted as 1 and controls as 0.}
  \item{data}{object of class "gwaa.data"}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
  \item{times}{number of times to do resampling analysis, 100 by default}
  \item{bcast}{progress is reported every \code{bcast} resumpling analyses}
}
\details{
In the analysis of empirical significance, first time the function 
\code{\link{ccfast}} is called and result object is 
saved. Later, the function \code{\link{ccfast}} is called 
\code{times} times with \code{resample=TRUE} in order to generate 
distribution under the null. Each call, minimal P-value is extracted 
and compared with original P-values. For a particular SNP, empirical 
P-value is obtained as a proportion of times minimal P from resampled data 
was less then original P.

The list elements effB, effAB and effBB are the ones obtained from the 
analysis of the original (not permuted) data set
}
\value{
	Object of class \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
\note{
The list elements medChi1df and medChi2df returned contain the median 
of the null distribution. 
}
\seealso{
\code{\link{ccfast}},
\code{\link{emp.qtscore}},
\code{\link{scan.gwaa-class}}
}
\examples{
data(srdta)
a<-ccfast("bt",data=srdta,snps=c(719:733))
plot(a)
# this does not make sense, as the whole experiment must be analysed, not a small region!
b<-emp.ccfast("bt",data=srdta,snps=c(719:733),bcast=10)
plot(b)
# compare qvalues and empirical P
qv<-qvaluebh95(a$P1df)$qval
qv
b$P1df
plot(qv,b$P1df)
abline(a=0,b=1)
}
\keyword{htest}% at least one, from doc/KEYWORDS
