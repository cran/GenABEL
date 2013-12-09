\name{polygenic_hglm}
\alias{polygenic_hglm}
\title{Estimation of polygenic model}
\usage{
  polygenic_hglm(formula, kinship.matrix, data,
    family = gaussian(), conv = 1e-06, maxit = 100, ...)
}
\arguments{
  \item{formula}{Formula describing fixed effects to be
  used in analysis, e.g. y ~ a + b means that outcome (y)
  depends on two covariates, a and b. If no covariates used
  in analysis, skip the right-hand side of the equation.}

  \item{kinship.matrix}{Kinship matrix, as provided by e.g.
  ibs(,weight="freq"), or estimated outside of GenABEL from
  pedigree data.}

  \item{data}{An (optional) object of
  \code{\link{gwaa.data-class}} or a data frame with
  outcome and covariates.}

  \item{family}{a description of the error distribution and
  link function to be used in the mean part of the model
  (see \code{\link{family}} for details of family
  functions).}

  \item{conv}{'conv' parameter of \code{hglm} (stricter
  than default, for great precision, use 1e-8).}

  \item{maxit}{'maxit' parameter of \code{hglm} (stricter
  than default).}

  \item{...}{other parameters passed to \code{hglm} call}
}
\description{
  Estimation of polygenic model using a hierarchical
  generalized linear model (HGLM; Lee and Nelder 1996.
  \code{hglm} package; Ronnegard et al. 2010).
}
\details{
  The algorithm gives extended quasi-likelihood (EQL)
  estimates for restricted maximum likelihood (REML)
  (Ronnegard et al. 2010). Implemented is the
  inter-connected GLM interpretation of HGLM (Lee and
  Nelder 2001). Testing on fixed and random effects
  estimates are directly done using inter-connected GLMs.
  Testing on dispersion parameters (variance components)
  can be done by extracting the profile likelihood function
  value of REML.
}
\examples{
require(GenABEL.data)
data(ge03d2ex.clean)
set.seed(1)
df <- ge03d2ex.clean[,sample(autosomal(ge03d2ex.clean),1000)]
gkin <- ibs(df,w="freq")

# ----- for quantitative traits
h2ht <- polygenic_hglm(height ~ sex + age, kin = gkin, df)
# ----- estimate of heritability
h2ht$esth2
# ----- other parameters
h2ht$h2an

# ----- test the significance of fixed effects
# ----- (e.g. quantitative trait)
h2ht <- polygenic_hglm(height ~ sex + age, kin = gkin, df)
# ----- summary with standard errors and p-values
summary(h2ht$hglm)
# ----- output for the fixed effects part
# ...
# MEAN MODEL
# Summary of the fixed effects estimates
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) 169.53525    2.57624  65.807  < 2e-16 ***
# sex          14.10694    1.41163   9.993  4.8e-15 ***
# age          -0.15332    0.05208  -2.944  0.00441 **
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Note: P-values are based on 69 degrees of freedom
# ...
# ----- extract fixed effects estimates and standard errors
h2ht$h2an

# ----- test the significance of polygenic effect
nullht <- lm(height ~ sex + age, df)
l1 <- h2ht$ProfLogLik
l0 <- as.numeric(logLik(nullht))
# the likelihood ratio test (LRT) statistic
(S <- -2*(l0 - l1))
# 5\% right-tailed significance cutoff
# for 50:50 mixture distribution between point mass 0 and \\chi^{2}(1)
# Self, SG and Liang KY (1987) Journal of the American Statistical Association.
qchisq(((1 - .05) - .50)/.50, 1)
# p-value
pchisq(S, 1, lower.tail = FALSE)/2

# ----- for binary traits
h2dm <- polygenic_hglm(dm2 ~ sex + age, kin = gkin, df, family = binomial(link = 'logit'))
# ----- estimated parameters
h2dm$h2an
}
\author{
  Xia Shen, Yurii Aulchenko
}
\references{
  Ronnegard, L, Shen, X and Alam, M (2010). hglm: A Package
  For Fitting Hierarchical Generalized Linear Models.
  \emph{The R Journal}, \bold{2}(2), 20-28.

  Lee, Y and Nelder JA (2001). Hierarchical generalised
  linear models: A synthesis of generalised linear models,
  random-effect models and structured dispersions.
  \emph{Biometrika}, \bold{88}(4), 987-1006.

  Lee, Y and Nelder JA (1996). Hierarchical Generalized
  Linear Models. \emph{Journal of the Royal Statistical
  Society. Series B}, \bold{58}(4), 619-678.
}
\seealso{
  \code{\link{polygenic}}, \code{\link{mmscore}},
  \code{\link{grammar}}
}
\keyword{htest}

