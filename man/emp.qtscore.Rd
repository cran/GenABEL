\name{emp.qtscore}
\alias{emp.qtscore}
\title{Bootstrap analysis of case-control GWA scan}
\description{
Bootstrap analysis of case-control GWA scan. Analysis 
function is \code{\link{qtscore}}. 
}
\usage{
emp.qtscore(formula , data, strata, snpsubset, idsubset, times = 100, 
		bcast = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  Only last two arguments are specific for \code{emp.qtscore} and others are 
  passed intact to \code{\link{qtscore}}.

    \item{formula}{Formula to be used in analysis. Should be a character string 
	following standard notation. On the left-had side, there should be outcome. 
	On the right-hand side, covariates are liste, with "+" separating the 
	covariates (additive action). The left- and right-hand sides are separated by "~". 
	For example "qt3~" would analyse asociation between SNPs and trait "qt3", 
	without any adjustment. To adjust for age and sex, use "qt3~age+sex".
	}
  \item{data}{object of class "gwaa.data"}
  \item{strata}{Stratification variable. If provieded, scores are computed within strata and 
		then added up.}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{times}{number of times to do resampling analysis, 100 by default}
  \item{bcast}{progress is reported every \code{bcast} resumpling analyses}
}
\details{
	In the analysis of empirical significance, first time the function 
	\code{\link{qtscore}} is called and result object is 
	saved. Later, the function \code{\link{qtscore}} is called 
	\code{times} times with \code{resample=TRUE} in order to generate 
	distribution under the null. Each call, minimal P-value is extracted 
	and compared with original P-values. For a particular SNP, empirical 
	P-value is obtained as a proportion of times minimal P from resampled data 
	was less then original P.

	The list elements effB, effAB and effBB are the ones obtained from the 
	analysis of the original (not permuted) data set

	The function does not yet implement correct analysis for X-linked data.
}
\value{
	Object of class \code{\link{scan.gwaa-class}}
}
\references{}
\author{Yurii Aulchenko}
\note{
The list elements medChi1df and medChi2df returned contain the median 
of the null distribution. 
}
\seealso{
\code{\link{qtscore}},
\code{\link{emp.ccfast}},
\code{\link{scan.gwaa-class}}
}
\examples{
data(srdta)
a<-qtscore("qt3~age+sex+CRSNP",data=srdta,snps=c(1:100))
plot(a)
b<-emp.qtscore("qt3~age+sex+CRSNP",data=srdta,snps=c(1:100))
plot(b)
}
\keyword{htest}% at least one, from doc/KEYWORDS
