\name{qtscore}
\alias{qtscore}
\title{Fast score test for association}
\usage{
  qtscore(formula, data, snpsubset, idsubset, strata,
    trait.type = "gaussian", times = 1, quiet = FALSE,
    bcast = 10, clambda = TRUE, propPs = 1, details = TRUE)
}
\arguments{
  \item{formula}{Formula describing fixed effects to be
  used in analysis, e.g. y ~ a + b means that outcome (y)
  depends on two covariates, a and b. If no covariates used
  in analysis, skip the right-hand side of the equation.}

  \item{data}{An object of \code{\link{gwaa.data-class}}}

  \item{snpsubset}{ndex, character or logical vector with
  subset of SNPs to run analysis on. If missing, all SNPs
  from \code{data} are used for analysis.}

  \item{idsubset}{ndex, character or logical vector with
  subset of IDs to run analysis on. If missing, all people
  from \code{data/cc} are used for analysis.}

  \item{strata}{Stratification variable. If provieded,
  scores are computed within strata and then added up.}

  \item{trait.type}{"gaussian" or "binomial" or "guess"
  (later option guesses trait type)}

  \item{times}{If more than one, the number of replicas to
  be used in derivation of empirical genome-wide
  significance. See \code{\link{emp.qtscore}}, which calls
  qtscore with times>1 for details}

  \item{quiet}{do not print warning messages}

  \item{bcast}{If the argument times > 1, progress is
  reported once in bcast replicas}

  \item{clambda}{If inflation facot Lambda is estimated as
  lower then one, this parameter controls if the original
  P1df (clambda=TRUE) to be reported in Pc1df, or the
  original 1df statistics is to be multiplied onto this
  "deflation" factor (clambda=FALSE). If a numeric value is
  provided, it is used as a correction factor.}

  \item{propPs}{proportion of non-corrected P-values used
  to estimate the inflation factor Lambda, passed directly
  to the \code{\link{estlambda}}}

  \item{details}{when FALSE, SNP and ID names are not
  reported in the returned object (saves some memory). This
  is experimental and will be not mantained anymore as soon
  as we achieve better memory efficiency for storage of SNP
  and ID names (currently default R character data type
  used)}
}
\value{
  Object of class \code{\link{scan.gwaa-class}}
}
\description{
  Fast score test for association between a trait and
  genetic polymorphism
}
\details{
  When formula contains covariates, the traits is analysed
  using GLM and later residuals used when score test is
  computed for each of the SNPs in analysis. Coefficients
  of regression are reported for the quantitative trait.

  For binary traits, odds ratios (ORs) are reportted. When
  adjustemnt is performed, first, "response" residuals are
  estimated after adjustment for covariates and scaled to
  [0,1]. Reported effects are approximately equal to ORs
  expected in logistic regression model.

  With no adjustment for binary traits, 1 d.f., the test is
  equivalent to the Armitage test.

  This is a valid function to analyse GWA data, including X
  chromosome. For X chromosome, stratified analysis is
  performed (strata=sex).
}
\examples{
require(GenABEL.data)
data(srdta)
#qtscore with stratification
a <- qtscore(qt3~sex,data=srdta)
plot(a)
b <- qtscore(qt3,strata=phdata(srdta)$sex,data=srdta)
add.plot(b,col="green",cex=2)
# qtscore with extra adjustment
a <- qtscore(qt3~sex+age,data=srdta)
a
plot(a)
# compare results of score and chi-square test for binary trait
a1 <- ccfast("bt",data=srdta,snps=c(1:100))
a2 <- qtscore(bt,data=srdta,snps=c(1:100),trait.type="binomial")
plot(a1,ylim=c(0,2))
add.plot(a2,col="red",cex=1.5)
# the good thing about score test is that we can do adjustment...
a2 <- qtscore(bt~age+sex,data=srdta,snps=c(1:100),trait.type="binomial")
points(a2[,"Position"],-log10(a2[,"P1df"]),col="green")
}
\author{
  Yurii Aulchenko
}
\references{
  Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid
  association using mixed model and regression: a fast and
  simple method for genome-wide pedigree-based quantitative
  trait loci association analysis. Genetics. 2007
  177(1):577-85.

  Amin N, van Duijn CM, Aulchenko YS. A genomic background
  based method for association analysis in related
  individuals. PLoS ONE. 2007 Dec 5;2(12):e1274.
}
\seealso{
  \code{\link{mlreg}}, \code{\link{mmscore}},
  \code{\link{egscore}}, \code{\link{emp.qtscore}},
  \code{\link{plot.scan.gwaa}},
  \code{\link{scan.gwaa-class}}
}
\keyword{htest}

