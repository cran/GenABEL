\name{PGC}
\alias{PGC}
\title{Polynomial genomic control}
\usage{
  PGC(data, method = "group_regress", p, df, pol.d = 3,
    plot = TRUE, index.filter = NULL, start.corr = FALSE,
    proportion = 1, n_quiantile = 5, title_name = "Lambda",
    type_of_plot = "plot", lmax = NULL, color = "red")
}
\arguments{
  \item{data}{Input vector of Chi square statistic}

  \item{method}{Function of error to be optimized. Can be
  "regress", "median", "ks.test" or "group_regress"}

  \item{p}{Input vector of allele frequencies}

  \item{df}{Number of degrees of freedom}

  \item{pol.d}{The degree of polinomial function}

  \item{plot}{If TRUE, plot of lambda will be produced}

  \item{start.corr}{For regress method use it only when you
  want to make calculations faster}

  \item{index.filter}{Index of variables in data vector,
  that will be used in analysis if zero - all variables
  will be used}

  \item{proportion}{The proportion of lowest P (Chi2) to be
  used when estimating the inflation factor Lambda for
  "regress" method only}

  \item{n_quiantile}{The number of groups for
  "group_regress" method}

  \item{title_name}{The title name for plot}

  \item{type_of_plot}{For developers only}

  \item{lmax}{The threshold for lambda for plotting
  (optional)}

  \item{color}{The color of the plot}
}
\value{
  A list with elements \item{data}{Output vector corrected
  Chi square statistic} \item{b}{Polinomial coefficients}
}
\description{
  This function estimates the genomic controls for
  different models and degrees of freedom, using polinomial
  function. Polinomial coefficients are estimated by
  optimizing different error functions: regress, median,
  ks.test or group regress.
}
\examples{
require(GenABEL.data)
data(ge03d2)
ge03d2 <- ge03d2[seq(from=1,to=nids(ge03d2),by=2),seq(from=1,to=nsnps(ge03d2),by=3)]
qts <- mlreg(dm2~1,data=ge03d2,gtmode = "additive")
chi2.1df <- results(qts)$chi2.1df
s <- summary(ge03d2)
freq <- s$Q.2
result=PGC(data=chi2.1df,method="median",p=freq,df=1, pol.d=2, plot=TRUE, lmax=1.1,start.corr=FALSE)
}
\author{
  Yakov Tsepilov
}
\keyword{htest}

