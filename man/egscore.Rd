\name{egscore}
\alias{egscore}
\title{Fast score test for association, corrected with PC}
\usage{
  egscore(formula, data, snpsubset, idsubset,
    kinship.matrix, naxes = 3, strata, times = 1,
    quiet = FALSE, bcast = 10, clambda = TRUE, propPs = 1)
}
\arguments{
  \item{formula}{Formula describing fixed effects to be
  used in analysis, e.g. y ~ a + b means that outcome (y)
  depends on two covariates, a and b.  If no covariates
  used in analysis, skip the right-hand side of the
  equation.}

  \item{data}{An object of \code{\link{gwaa.data-class}}}

  \item{snpsubset}{Index, character or logical vector with
  subset of SNPs to run analysis on.  If missing, all SNPs
  from \code{data} are used for analysis.}

  \item{idsubset}{Index, character or logical vector with
  subset of IDs to run analysis on.  If missing, all people
  from \code{data/cc} are used for analysis.}

  \item{kinship.matrix}{kinship matrix, as returned by
  \code{\link{ibs}}, Use weight="freq" with
  \code{\link{ibs}} and do not forget to repalce the
  diagonal with Var returned by \code{\link{hom}}, as shown
  in example!}

  \item{naxes}{Number of axes of variation to be used in
  adjustment (should be much smaller than number of
  subjects)}

  \item{strata}{Stratification variable. If provieded,
  scores are computed within strata and then added up.}

  \item{times}{If more then one, the number of replicas to
  be used in derivation of empirical genome-wide
  significance.}

  \item{quiet}{do not print warning messages}

  \item{bcast}{If the argument times > 1, progress is
  reported once in bcast replicas}

  \item{clambda}{If inflation facot Lambda is estimated as
  lower then one, this parameter controls if the original
  P1df (clambda=TRUE) to be reported in Pc1df, or the
  original 1df statistics is to be multiplied onto this
  "deflation" factor (clambda=FALSE).  If a numeric value
  is provided, it is used as a correction factor.}

  \item{propPs}{proportion of non-corrected P-values used
  to estimate the inflation factor Lambda, passed directly
  to the \code{\link{estlambda}}}
}
\value{
  Object of class \code{\link{scan.gwaa-class}}
}
\description{
  Fast score test for association between a trait and
  genetic polymorphism, adjusted for possible
  stratification by principal components.
}
\details{
  The idea of this test is to use genomic kinship matrix to
  first, derive axes of genetic variation (principal
  components), and, second, adjust both trait and genotypes
  onto these axes. Note that the diagonal of the kinship
  matrix should be replaced (default it is 0.5*(1+F), and
  for EIGENSTRAT one needs variance).  These variances are
  porduced by \code{\link{hom}} function (see example).

  The traits is first analysed using LM and with covariates
  as specified with formula and also with axes of variation
  as predictors. Corrected genotypes are defined as
  residuals from regression of genotypes onto axes (which
  are orthogonal). Correlaton between corrected genotypes
  and phenotype is computed, and test statistics is defined
  as square of this correlation times (N - K - 1), where N
  is number of genotyped subjects and K is the number of
  axes.

  This test is defined only for 1 d.f.
}
\examples{
require(GenABEL.data)
data(ge03d2c)
#egscore with stratification
gkin <- ibs(ge03d2c[,autosomal(ge03d2c)],w="freq")
#replace the diagonal with right elements
diag(gkin) <- hom(ge03d2c[,autosomal(ge03d2c)])$Var
a <- egscore(dm2~sex+age,data=ge03d2c,kin=gkin)
plot(a,df="Pc1df")
}
\author{
  Yurii Aulchenko
}
\references{
  Price A. L. et al, Principal components analysis corrects
  for stratification in genome-wide association studies.
  Nat Genet 38: 904-909.
}
\seealso{
  \code{\link{qtscore}}, \code{\link{mmscore}},
  \code{\link{ibs}}, \code{\link{scan.gwaa-class}}
}
\keyword{htest}

