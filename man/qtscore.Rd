\name{qtscore}
\alias{qtscore}
\title{fast score test for association}
\description{
Fast score test for association between a trait and genetic polymorphism
}
\usage{
qtscore(formula, data, strata, snpsubset, idsubset, quiet = FALSE) 
}
\arguments{
  \item{formula}{Formula to be used in analysis. It should be a character string 
	following standard notation. On the left-had side, there should be outcome. 
	On the right-hand side, covariates are liste, with "+" separating the 
	covariates (additive action). The left- and right-hand sides are separated by "~". 
	You should put CRSNP argument in the formula.
	For example "qt3~CRSNP" would analyse asociation between SNPs and trait "qt3", 
	without any adjustment. To adjust for age and sex, use "qt3~age+sex+CRSNP".
	At current stage, only additive effects ("+") are allowed.
	}
  \item{data}{object of class "gwaa.data"}
  \item{strata}{Stratification variable. If provieded, scores are computed within strata and 
		then added up.}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{quiet}{do not print warning messages}
}
\details{
When formula contains covariates, the traits is analysed using GLM and later 
residuals used when score test is computed for each of the SNPs in analysis. 
It appears that for simple stratas (e.g. sex), this adjustment and analysis by 
strata generate very similar results (see example).

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
add.plot.scan.gwaa(b,col="green",cex=2)
#atscore with extra adjustment
a <- qtscore("qt3~sex+age+CRSNP",data=srdta,snps=c(1:10))
a
plot(a)
# compare results of score and chi-square test for binary trait
a1 <- ccfast("bt",data=srdta,snps=c(1:100))
a2 <- qtscore("bt~CRSNP",data=srdta,snps=c(1:100))
plot(a1,ylim=c(0,2))
points(a2$map,-log10(a2$P1df),col="red")
# the good thing about score test is that we can do adjustment...
a2 <- qtscore("bt~age+sex+CRSNP",data=srdta,snps=c(1:100))
points(a2$map,-log10(a2$P1df),col="green")
}
\keyword{htest}% at least one, from doc/KEYWORDS
