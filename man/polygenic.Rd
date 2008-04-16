\name{polygenic}
\alias{polygenic}
\title{Estimation of polygenic model}
\description{
	Estimates linear mixed (polygenic) model based on trait and covariates data 
	and kinship matrix
}
\usage{
polygenic(formula,kinship.matrix,data,fixh2,starth2=0.3,trait.type="gaussian",opt.method="nlm",scaleh2=1000,quiet=FALSE,...)
}
\arguments{
  \item{formula}{Formula describing fixed effects to be used in analysis, e.g. 
	y ~ a + b means that outcome (y) depends on two covariates, a and b. 
	If no covariates used in analysis, skip the right-hand side of the 
	equation.
	}
  \item{kinship.matrix}{Kinship matrix, as provided by e.g. ibs(,weight="freq"), 
	or estimated outside of GenABEL from pedigree data.}
  \item{data}{An (optional) object of \code{\link{gwaa.data-class}} or a data frame with 
	outcome and covariates
	}
  \item{fixh2}{Optional value of heritability to be used, instead of maximisation. 
	The uses of this option are two-fold: (a) testing significance of 
	heritability and (b) using a priori known heritability to derive the 
	rest of MLEs and var.-cov. matrix.
	}
  \item{starth2}{Starting value for h2 estimate}
  \item{trait.type}{"gaussian" or "binomial"}
  \item{opt.method}{"nlm" or "optim". These two use dirrerent optimisation functions. 
	\code{\link{optim}} is slower than \code{\link{nlm}}, but may give better 
	results.
	}
  \item{scaleh2}{Only relevant when "nlm" optimisation function is used. 
	"scaleh2" is the heritability 
	scaling parameter, regulating how "big" are parameter changes in h2 with the 
	respect to changes in other parameters. As other parameters are estimated 
	from previous regression, these are expected to change little from the 
	initial estimate. The default value of 1000 proved to work rather well under a 
	range of conditions.
	}
  \item{quiet}{If FALSE (default), details of optimisation process are reported. 
	}
  \item{...}{Optional arguments to be passed to \code{\link{nlm}} (\code{\link{optim}}) 
	minimisation function.
	}
}
\details{
	This function maximises the likelihood of the data under polygenic 
	model with covariates an reports the maximum likelihood estimates 
	and the inverse of variance-covariance matrix at the point of ML. 

	One of the major use of this function is to estimate residuals of the 
	trait and the inverse of the variance-covariance matrix for 
	further use in analysis with \code{\link{mmscore}} and 
	\code{\link{grammar}}.

	Also, it can be used for a variant of GRAMMAR analysis, which 
	allows for permutations for GW significance by use of 
	environmental residuals as an analysis trait with \code{\link{qtscore}}.

	"Environmental residuals" (not to be mistaken with just "residuals") are 
	the residual where both the effect of covariates AND the estimated 
	polygenic effect (breeding values) are factored out. This thus 
	provides an estimate of the trait value contributed by environment
	(or, turning this other way around, the part of trait not explained 
	by covariates and by the polygene). Polygenic residuals are estimated 
	as

	$$
	\sigma^2  V^{-1} (Y - (\hat{\mu} + \hat{\beta} C_1 + ...))
	$$

	where $sigma^2$ is the residual variance, $V^{-1}$ is the 
	InvSigma (inverse of the var-cov matrix at the maximum of 
	polygenic model) and 
	$(Y - (\hat{\mu} + \hat{\beta} C_1 + ...))$ is the trait 
	values adjusted for covariates (also at at the maximum of 
	polygenic model likelihood). 

	It can also be used for heritability analysis.
	If you want to test significance of heritability, 
	estimate the model and write down 
	the function minimum reported at "h2an" element of the output 
	(this is -2*MaxLikleihood). Then do next round of 
	estimation, but set fixh2=0. The difference between you function minima 
	gives you one-sided test distribued as chi-squared with 1 d.f.

	The way to compute the likleihood is partly based on 
	the paper of Thompson (see refs), namely instead of 
	taking inverse of var-cov matrix every time, 
	eigenvectors of the inverse of G (taken only once) 
	are used.
}
\value{
  A list with values 
  \item{h2an}{A list supplied by the \code{\link{nlm}} minimisation routine. 
	Of particular interest are elements "estimate" containing parameter 
	maximal likelihood estimates (MLEs) (order: mean, betas for covariates, 
	heritability, (polygenic + residual variance))}
  \item{residualY}{Residuals from analysis, based on covariate effects only; 
	NOTE: these are NOT grammar "environmental residuals"!}
  \item{esth2}{Estimate (or fixed value) of heritability}
  \item{pgresidualY}{Environmental residuals from analysis, based on covariate effects 
	and predicted breeding value.
	}
  \item{InvSigma}{Inverse of the variance-covariance matrix, computed at the 
	MLEs -- these are used in \code{\link{mmscore}} and \code{\link{grammar}}
	functions.}
  \item{call}{The details of call}
  \item{measuredIDs}{Logical values for IDs who were used in analysis 
	(traits and all covariates measured) == TRUE}
}
\references{
Thompson EA, Shaw RG (1990) Pedigree analysis for quantitative 
traits: variance components without matrix inversion. Biometrics 
46, 399-413.

Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid association using mixed model 
and regression: a fast and simple method for genome-wide pedigree-based quantitative 
trait loci association analysis. Genetics. 2007 177(1):577-85.

Amin N, van Duijn CM, Aulchenko YS. A genomic background based method for 
association analysis in related individuals. PLoS ONE. 2007 Dec 5;2(12):e1274. 
}
\author{Yurii Aulchenko}
\note{
	Presence of twins may screw up your analysis. Check kinship matrix for 
	singularities, or rather use \code{\link{check.marker}} for identification 
	and exclusion of twin samples.

	If a trait (no covarites) is used, make sure that order of IDs in 
	kinship.matrix is exactly the same as in the outcome
}
\seealso{
\code{\link{mmscore}},
\code{\link{grammar}}
}
\examples{
# note that procedure runs on CLEAN data
data(ge03d2ex.clean)
gkin <- ibs(ge03d2ex.clean,w="freq")
h2ht <- polygenic(height ~ sex + age,kin=gkin,ge03d2ex.clean)
# estimate of heritability
h2ht$esth2
# other parameters
h2ht$h2an
# minimum negative log-likelihood
h2ht$h2an$minimum
# maximum log-likelihood
-h2ht$h2an$minimum

#for binary trait (experimental)
h2dm <- polygenic(dm2 ~ sex + age,kin=gkin,ge03d2ex.clean,trait="binomial")
# estimated parameters
h2dm$h2an
}
\keyword{htest}% at least one, from doc/KEYWORDS
