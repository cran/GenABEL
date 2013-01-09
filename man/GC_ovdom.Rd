\name{GC_ovdom}
\alias{GC_ovdom}
\title{Genomic control for overdimunant model of inheritance using VIF}
\usage{
  GC_ovdom(data, p, method = "regress", n,
    index.filter = NULL, proportion = 1, clust = 0,
    vart0 = 0, tmp = 0)
}
\arguments{
  \item{data}{Input vector of Chi square statistic}

  \item{method}{Function of error to be optimized. Can be
  "regress", "median" or "ks.test"}

  \item{p}{Input vector of allele frequencies}

  \item{index.filter}{Indexes for variables that will be
  use for analisis in data vector}

  \item{n}{size of the sample}

  \item{proportion}{The proportion of lowest P (Chi2) to be
  used when estimating the inflation factor Lambda for
  "regress" method only}

  \item{clust}{For developers only}

  \item{vart0}{For developers only}

  \item{tmp}{For developers only}
}
\value{
  A list with elements \item{Zx}{output vector corrected
  Chi square statistic} \item{vv}{output vector of VIF}
  \item{exeps}{output vector of exepsons (NA)}
  \item{calrate}{output vector of calrate} \item{F}{F}
  \item{K}{K}
}
\description{
  This function estimates the corrected statistic using
  genomic control for overdimunant model, using VIF. VIF
  coefficients are estimated by optimizing diffrent error
  functions: regress, median and ks.test.
}
\examples{
data(ge03d2)
# truncate the data to make the example faster
ge03d2 <- ge03d2[seq(from=1,to=nids(ge03d2),by=2),seq(from=1,to=nsnps(ge03d2),by=3)]
qts <- mlreg(phdata(ge03d2)$dm2~1,data=ge03d2,gtmode = "overdominant")
chi2.1df <- results(qts)$chi2.1df
s <- summary(ge03d2)
freq <- s$Q.2
result <- GC_ovdom(p=freq,method = "median",data=chi2.1df,n=nids(ge03d2))
}
\author{
  Yakov Tsepilov
}
\keyword{htest}

