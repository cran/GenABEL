\name{PGC}
\alias{PGC}
\title{Polynomial genomic control}
\usage{
  PGC(data, method = "regress", p, df, pol.d = 3,
    plot = TRUE, index.filter = NULL, start.corr = FALSE,
    proportion = 1)
}
\arguments{
  \item{data}{Input vector of Chi square statistic}

  \item{method}{Function of error to be optimized. Can be
  "regress", "median" or "ks.test"}

  \item{p}{Input vector of allele frequencies}

  \item{df}{Number of dergees of freedom}

  \item{pol.d}{Polinom degree}

  \item{plot}{If true, functon makes plot of lambda from
  allele frequencies}

  \item{start.corr}{For regress method use it only when you
  want to make calculations faster}

  \item{index.filter}{Index of variables in data vector,
  that will be used in analysis if zero - all variables
  will be used}

  \item{proportion}{The proportion of lowest P (Chi2) to be
  used when estimating the inflation factor Lambda for
  "regress" method only}
}
\value{
  A list with elements \item{data}{Output vector corrected
  Chi square statistic} \item{b}{Polinom coefficients}
}
\description{
  This function estimates the genomic controls for diffrent
  model and degrees of freedom, using polinom function.
  Polinomic coefficients are estimated by optimizing
  diffrent error functions: regress, median and ks.test.
}
\examples{
data(ge03d2)
qts=qtscore(dm2, ge03d2)
chi2.1df=results(qts)$chi2.1df
s=summary(ge03d2)
MAF=s$Q.2
result=PGC(data=chi2.1df,method="regress",p=MAF,df=1, pol.d=3, plot=FALSE, start.corr=FALSE)
}
\author{
  Yakov Tsepilov
}
\keyword{htest}

