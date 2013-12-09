\name{mlreg.p}
\alias{mlreg.p}
\title{EXPERIMENTAL Linear and logistic regression and Cox models for genome-wide SNP data}
\description{
Linear and logistic regression and Cox models for genome-wide SNP data
}
\usage{
mlreg.p(formula, data, snpsubset, idsubset, gtmode = "additive", trait.type = "guess")
%mlreg.p(formula, data, gtmode = "additive", trait.type = "guess", propPs = 1)
}
\arguments{
  \item{formula}{Standard formula object}
  \item{data}{an object of \code{\link{gwaa.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{gtmode}{Either "additive", "dominant", "recessive" or "overdominant". Specifies the 
		analysis model.}
  \item{trait.type}{Either "gaussian", "binomial" or "survival", corresponding to 
		analysis using linear regression, logistic regression, and Cox proportional 
		hazards models, respectively. When default vale "guess" is used, the 
		program tries to guess the type}
}
\details{
	Linear regression is performed using standard approach; logisitc regression 
	is implemented using IRLS; Cox model makes use of code contributed by 
	Thomas Lumley (survival package).

	For logistic and Cox, exp(effB) gives Odds Ratios and Hazard Ratios, respectively.
}
\value{
	An object of \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{ 
\code{\link{GASurv}},
\code{\link{qtscore}}
}
\examples{
	require(GenABEL.data)
	data(ge03d2)
	dta <- ge03d2[,1:100]
# analysis using linear model
	xq <- mlreg.p(bmi~sex,dta)
# logistic regression, type guessed automatically
	xb <- mlreg.p(dm2~sex,dta)
# Cox proportional hazards model, assuming that age is the follow-up time 
# generally this does not make sense (could be ok if age is age at onset)
	xs <- mlreg.p(GASurv(age,dm2)~sex,dta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
