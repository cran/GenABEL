\name{qtscore}
\alias{qtscore}
\title{Fast score test for association}
\description{
Fast score test for association between a trait and genetic polymorphism
}
\usage{
qtscore(formula,data,snpsubset,idsubset,strata,trait.type,times=1,quiet=FALSE,bcast=10,clambda=TRUE)
}
\arguments{
  \item{formula}{Formula to be used in analysis. It should be a character string 
	(thus quotation marks from both sides) 
	following the standard notation. On the left-had side, there should be outcome. 
	On the right-hand side, covariates are liste, with "+" separating the 
	covariates (additive action). The left- and right-hand sides are separated by "~". 
	You should put CRSNP argument in the formula.
	For example "qt3~CRSNP" assumes a crude analysis of association between SNPs and 
	trait "qt3". To adjust for e.g. age and sex, use "qt3~age+sex+CRSNP".
	At current stage, only additive effects ("+") are allowed.
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
  \item{clambda}{If inflation facot Lambda is estimated as lower then one, this parameter 
		controls if the original P1df (clambda=TRUE) to be reported in Pc1df, 
		or the original 1df statistics is to be multiplied onto this "deflation" 
		facor (clambda=FALSE)}
}
\details{
When formula contains covariates, the traits is analysed using GLM and later 
residuals used when score test is computed for each of the SNPs in analysis. 
For binary traits, residuals from GLM are transformed using exp(x)/(1+exp(x)).

With no adjustment for binary traits, 1 d.f., the test is equivalent to the 
Armitage test.

This is a valid function to analyse GWA data, including X chromosome. For X chromosome,
stratified analysis is performed (strata=sex).
}
\value{
  Object of class \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{emp.qtscore}},
\code{\link{plot.scan.gwaa}},
\code{\link{scan.gwaa-class}}
}
\examples{
data(srdta)
#qtscore with stratification
a <- qtscore("qt3~sex+CRSNP",data=srdta)
plot(a)
b <- qtscore("qt3~CRSNP",strata=srdta@phdata$sex,data=srdta)
add.plot(b,col="green",cex=2)
# qtscore with extra adjustment
a <- qtscore("qt3~sex+age+CRSNP",data=srdta)
a
plot(a)
# compare results of score and chi-square test for binary trait
a1 <- ccfast("bt",data=srdta,snps=c(1:100))
a2 <- qtscore("bt~CRSNP",data=srdta,snps=c(1:100))
plot(a1,ylim=c(0,2))
add.plot(a2,col="red",cex=1.5)
# the good thing about score test is that we can do adjustment...
a2 <- qtscore("bt~age+sex+CRSNP",data=srdta,snps=c(1:100))
points(a2$map,-log10(a2$P1df),col="green")
}
\keyword{htest}% at least one, from doc/KEYWORDS
