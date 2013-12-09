\name{polygenic}
\alias{polygenic}
\title{Estimation of polygenic model}
\usage{
  polygenic(formula, kinship.matrix, data, fixh2,
    starth2 = 0.3, trait.type = "gaussian",
    opt.method = "nlm", scaleh2 = 1, quiet = FALSE,
    steptol = 1e-08, gradtol = 1e-08, optimbou = 8,
    fglschecks = TRUE, maxnfgls = 8, maxdiffgls = 1e-04,
    patchBasedOnFGLS = TRUE, llfun = "polylik_eigen",
    eigenOfRel, ...)
}
\arguments{
  \item{formula}{Formula describing fixed effects to be
  used in the analysis, e.g. y ~ a + b means that outcome
  (y) depends on two covariates, a and b. If no covariates
  used in the analysis, skip the right-hand side of the
  equation.}

  \item{kinship.matrix}{Kinship matrix, as provided by e.g.
  ibs(,weight="freq"), or estimated outside of GenABEL from
  pedigree data.}

  \item{data}{An (optional) object of
  \code{\link{gwaa.data-class}} or a data frame with
  outcome and covariates}

  \item{fixh2}{Optional value of heritability to be used,
  instead of maximisation. The uses of this option are
  two-fold: (a) testing significance of heritability and
  (b) using a priori known heritability to derive the rest
  of MLEs and var.-cov. matrix.}

  \item{starth2}{Starting value for h2 estimate}

  \item{trait.type}{"gaussian" or "binomial"}

  \item{opt.method}{"nlm" or "optim". These two use
  different optimisation functions. We suggest using the
  default \code{\link{nlm}}, although \code{\link{optim}}
  may give better results in some situations}

  \item{scaleh2}{Only relevant when "nlm" optimisation
  function is used. "scaleh2" is the heritability scaling
  parameter, regulating how "big" are parameter changes in
  h2 with respect to changes in other parameters. As other
  parameters are estimated from previous regression, these
  are expected to change little from the initial estimate.
  The default value of 1000 proved to work rather well
  under a range of conditions.}

  \item{quiet}{If FALSE (default), details of optimisation
  process are reported}

  \item{steptol}{steptal parameter of "nlm"}

  \item{gradtol}{gradtol parameter of "nlm"}

  \item{optimbou}{fixed effects boundary scale parameter
  for 'optim'}

  \item{fglschecks}{additional check for convergence on/off
  (convergence between estimates obtained and that from
  FGLS)}

  \item{maxnfgls}{number of fgls checks to perform}

  \item{maxdiffgls}{max difference allowed in fgls checks}

  \item{patchBasedOnFGLS}{if FGLS checks not passed,
  'patch' fixed effect estimates based on FGLS expectation}

  \item{llfun}{function to compute likelihood (default
  'polylik_eigen', also available -- but not recommended --
  'polylik')}

  \item{eigenOfRel}{results of eigen(relationship matrix =
  2*kinship.matrix). Passing this can decrease
  computational time substantially if multiple traits are
  analysed using the same kinship matrix. This option will
  not work if any NA's are found in the trait and/or
  covariates and if the dimensions of the 'eigen'-object,
  trait, covariates, kinship do not match.}

  \item{...}{Optional arguments to be passed to
  \code{\link{nlm}} or (\code{\link{optim}}) minimisation
  function}
}
\value{
  A list with values \item{h2an}{A list supplied by the
  \code{\link{nlm}} minimisation routine. Of particular
  interest are elements "estimate" containing parameter
  maximal likelihood estimates (MLEs) (order: mean, betas
  for covariates, heritability, (polygenic + residual
  variance)). The value of twice negative maximum
  log-likelihood is returned as h2an\$minimum.}
  \item{esth2}{Estimate (or fixed value) of heritability}
  \item{residualY}{Residuals from analysis, based on
  covariate effects only; NOTE: these are NOT grammar
  "environmental residuals"!}
  \item{pgresidualY}{Environmental residuals from analysis,
  based on covariate effects and predicted breeding value.
  } \item{grresidualY}{GRAMMAR+ transformed trait
  residuals} \item{grammarGamma}{list with GRAMMAR-gamma
  correction factors} \item{InvSigma}{Inverse of the
  variance-covariance matrix, computed at the MLEs -- these
  are used in \code{\link{mmscore}} and
  \code{\link{grammar}} functions.} \item{call}{The details
  of call} \item{measuredIDs}{Logical values for IDs who
  were used in analysis (traits and all covariates
  measured) == TRUE} \item{convFGLS}{was convergence
  achieved according to FGLS criterionE}
}
\description{
  This function maximises the likelihood of the data under
  polygenic model with covariates an reports twice negative
  maximum likelihood estimates and the inverse of the
  variance-covariance matrix at the point of ML.
}
\details{
  One of the major uses of this function is to estimate
  residuals of the trait and the inverse of the
  variance-covariance matrix for further use in analysis
  with \code{\link{mmscore}} and \code{\link{grammar}}.

  Also, it can be used for a variant of GRAMMAR analysis,
  which allows for permutations for GW significance by use
  of environmental residuals as an analysis trait with
  \code{\link{qtscore}}.

  "Environmental residuals" (not to be mistaken with just
  "residuals") are the residual where both the effect of
  covariates AND the estimated polygenic effect (breeding
  values) are factored out. This thus provides an estimate
  of the trait value contributed by environment (or,
  turning this other way around, the part of the trait not
  explained by covariates and by the polygene). Polygenic
  residuals are estimated as

  \deqn{ \sigma^2 V^{-1} (Y - (\hat{\mu} + \hat{\beta} C_1
  + ...)) }

  where \eqn{sigma^2} is the residual variance,
  \eqn{V^{-1}} is the InvSigma (inverse of the var-cov
  matrix at the maximum of polygenic model) and \eqn{(Y -
  (\hat{\mu} + \hat{\beta} C_1 + ...))} is the trait values
  adjusted for covariates (also at at the maximum of
  polygenic model likelihood).

  It can also be used for heritability analysis. If you
  want to test significance of heritability, estimate the
  model and write down the function minimum reported at the
  "h2an" element of the output (this is twice the negative
  MaxLikelihood). Then do a next round of estimation, but
  set fixh2=0. The difference between your function minima
  gives a test distributed as chi-squared with 1 d.f.

  The way to compute the likelihood is partly based on the
  paper of Thompson (see refs), namely instead of taking
  the inverse of the var-cov matrix every time,
  eigenvectors of the inverse of G (taken only once) are
  used.
}
\note{
  Presence of twins may complicate your analysis. Check the
  kinship matrix for singularities, or rather use
  \code{\link{check.marker}} for identification of twin
  samples. Take special care in interpretation.

  If a trait (no covariates) is used, make sure that the
  order of IDs in the kinship.matrix is exactly the same as
  in the outcome

  Please note that there is alternative to 'polygenic',
  \code{\link{polygenic_hglm}}, which is faster than
  polygenic() with the llfun='polylik' option, but slightly
  slower than the default polygenic().
}
\examples{
# note that procedure runs on CLEAN data
require(GenABEL.data)
data(ge03d2ex.clean)
gkin <- ibs(ge03d2ex.clean,w="freq")
h2ht <- polygenic(height ~ sex + age, kin=gkin, ge03d2ex.clean)
# estimate of heritability
h2ht$esth2
# other parameters
h2ht$h2an
# the minimum twice negative log-likelihood
h2ht$h2an$minimum
# twice maximum log-likelihood
-h2ht$h2an$minimum

# for binary trait (experimental)
h2dm <- polygenic(dm2 ~ sex + age, kin=gkin, ge03d2ex.clean, trait="binomial")
# estimated parameters
h2dm$h2an
}
\author{
  Yurii Aulchenko, Gulnara Svischeva
}
\references{
  Thompson EA, Shaw RG (1990) Pedigree analysis for
  quantitative traits: variance components without matrix
  inversion. Biometrics 46, 399-413.

  for original GRAMMAR

  Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid
  association using mixed model and regression: a fast and
  simple method for genome-wide pedigree-based quantitative
  trait loci association analysis. Genetics. 2007
  177(1):577-85.

  for GRAMMAR-GC

  Amin N, van Duijn CM, Aulchenko YS. A genomic background
  based method for association analysis in related
  individuals. PLoS ONE. 2007 Dec 5;2(12):e1274.

  for GRAMMAR-Gamma

  Svischeva G, Axenovich TI, Belonogova NM, van Duijn CM,
  Aulchenko YS. Rapid variance components-based method for
  whole-genome association analysis. Nature Genetics. 2012
  44:1166-1170. doi:10.1038/ng.2410

  for GRAMMAR+ transformation

  Belonogova NM, Svishcheva GR, van Duijn CM, Aulchenko YS,
  Axenovich TI (2013) Region-Based Association Analysis of
  Human Quantitative Traits in Related Individuals. PLoS
  ONE 8(6): e65395. doi:10.1371/journal.pone.0065395
}
\seealso{
  \code{\link{polygenic_hglm}}, \code{\link{mmscore}},
  \code{\link{grammar}}
}
\keyword{htest}

