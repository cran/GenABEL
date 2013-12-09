\name{emp.qtscore}
\alias{emp.qtscore}
\title{Genome-wide significance for a GWA scan}
\description{
Genome-wide significance for a GWA scan. Analysis 
function is \code{\link{qtscore}}. 
}
\usage{
emp.qtscore(formula , data, snpsubset, idsubset, strata, trait.type="gaussian", 
		times = 200, quiet=FALSE, bcast = 10)
}
\arguments{
  All arguments are the same as in and passed intact to the \code{\link{qtscore}}.
  See help for this function. 
  \item{formula}{Formula describing fixed effects to be used in analysis, e.g. 
	y ~ a + b means that outcome (y) depends on two covariates, a and b. 
	If no covariates used in analysis, skip the right-hand side of the 
	equation.
	}
  \item{data}{An object of \code{\link{gwaa.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{strata}{Stratification variable. If provieded, scores are computed within strata and 
		then added up.}
  \item{trait.type}{"gaussian" or "binomial". If not specified, the procedure quesses the type}
  \item{times}{If more then one, the number of replicas to be used in derivation of 
		empirical genome-wide significance. See \code{\link{emp.qtscore}}, which
		calls qtscore with times>1 for details}
  \item{quiet}{do not print warning messages}
  \item{bcast}{If the argument times > 1, progress is reported once in bcast replicas}
}
\details{
	In the analysis of empirical significance, first time the function 
	\code{\link{qtscore}} is called and result object is 
	saved. Later, the function \code{\link{qtscore}} is called 
	\code{times} times with \code{replace=FALSE} in order to generate 
	distribution under the null. Each call, minimal P-value is extracted 
	and compared with original P-values. For a particular SNP, empirical 
	P-value is obtained as a proportion of times minimal Ps from resampled data 
	was less then original P.

	The list elements effB, effAB and effBB are the ones obtained from the 
	analysis of the original (not permuted) data set

	The function does not yet implement correct analysis for X-linked data.
}
\value{
	Object of class \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{qtscore}},
\code{\link{emp.ccfast}},
\code{\link{scan.gwaa-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
a<-qtscore(qt3~age+sex,data=srdta,snps=c(1:200))
plot(a)
# this does not make sense, as the whole experiment must be analysed, not a small region!
# also, times = 20 is way too small (should be at least 200)
b<-emp.qtscore(qt3~age+sex,data=srdta,snps=c(1:200), times = 20)
plot(b)
}
\keyword{htest}% at least one, from doc/KEYWORDS
