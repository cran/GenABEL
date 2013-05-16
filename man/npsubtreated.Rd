\name{npsubtreated}
\alias{npsubtreated}
\title{non-parametric trait "imputations" in treated people}
\description{
	For people on treatment, the algorithm substitutes the value of the trait 
	using non-parametric algorithm described in Tobin et al., 2005. This algorithm 
	assumes that the measurement in treated subject is a right-censored trait. 
	Essentially, the algorithm substitutes the QT for a person on treatment with 
	the mean of the above-ranked substituted QT value.
}
\usage{
npsubtreated(trait, medication, increase = FALSE)
}
\arguments{
  \item{trait}{vector if trait values}
  \item{medication}{medication indicator (0/1)}
  \item{increase}{if medication INCREASE the value of the trait (should never be true 
	for e.g. blood pressure, LDL cholesterol, etc.)}
}
\details{
	Should put the formulas here...
}
\value{
	Vector of trait walues, where the values for treated people are substituted 
	with average of the above-ranked substituted trait value.
}
\references{
Levy D, DeStefano AL, Larson MG, O'Donnell CJ, Lifton RP, Gavras H, Cupples LA, 
Myers RH. Evidence for a gene influencing blood pressure on chromosome 17. 
Genome scan linkage results for longitudinal blood pressure phenotypes in 
subjects from the framingham heart study. Hypertension. 2000 Oct;36(4):477-83.

Tobin MD, Sheehan NA, Scurrah KJ, Burton PR.
Adjusting for treatment effects in studies of quantitative traits: 
antihypertensive therapy and systolic blood pressure.
Stat Med. 2005 Oct 15;24(19):2911-35. 
}
\author{Yurii Aulchenko}
%\note{
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
# simulate SBP data
	simmeddat <- function(mu=144,bage=0.5,bsex=4.,bg=2.,pB=0.3,rvar=21^2,N=1000) {
		ageb <- c(25,74)
		pmale <- .5
		htthresh <- 160
		trprob <- .5
		mutreff <- (-15.)
		trvar <- 4^2
		age <- runif(N,min=ageb[1],max=ageb[2])
		sex <- 1*(runif(N)<=pmale)
		gt <- rbinom(N,size=2,prob=pB)
		y.true <- rnorm(N,mu,sqrt(rvar)) + bage*age + bsex*sex + bg*gt
		d.true <- (y.true>=htthresh)
		medication <- 1*d.true
		medication[d.true] <- 1*(runif(sum(d.true))<=trprob)
		treatm <- rnorm(sum(medication),mutreff,sqrt(trvar))
		treatm[treatm>0] <- 0
		treff <- rep(0,N)
		treff[medication==1] <- treatm
		y.obs <- y.true + treff
		out <- data.frame(age,sex,gt,y.true,d.true,medication,treff,y.obs)
		out
}
x <- simmeddat(bg=2.0,N=3000)
x[1:15,]

# substitute value of treated people
imptra <- npsubtreated(x$y.obs,x$medication)
imptra[1:15]

# Almost always, correlation should be higher for the "imputed" trait
cor(x$y.true,x$y.obs)
cor(x$y.true,imptra)

# see what comes out of regression
# analysis of true value
summary(lm(y.true~sex+age+gt,data=x))
# ignore treatment (as a rule, betas are underestimated; 
# on average, power goes down)
summary(lm(y.obs~sex+age+gt,data=x))
# treatment as covariate (as a rule, betas are underestimated; 
# on average, power goes down)
summary(lm(y.obs~sex+age+gt+medication,data=x))
# analyse "imputed" trait (as a rule betas are better; power 
# approaches that of analysis of "true" trait)
summary(lm(imptra~sex+age+gt,data=x))
}
\keyword{robust}
