\name{scan.haplo.2D}
\alias{scan.haplo.2D}
\title{runs haplo.score.slide with all pairs of markers in a region}
\description{
	Runs \code{\link[haplo.stats]{haplo.score.slide}} from the package 
	\code{haplo.stats} on all pairs of markers in a region 
	and presents output as \code{\link{scan.gwaa.2D-class}} object
}
\usage{
scan.haplo.2D(formula, data, snpsubset, idsubset, bcast = 10, 
				simulate=FALSE, trait.type, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Formula to be used in analysis. It should be a character string 
	following standard notation. On the left-had side, there should be outcome. 
	On the right-hand side, covariates are liste, with "+" separating the 
	covariates (additive action). The left- and right-hand sides are separated by "~". 
	You should put CRSNP argument in the formula.
	For example "qt3~CRSNP" would analyse asociation between SNPs and trait "qt3", 
	without any adjustment. To adjust for age and sex, use "qt3~age+sex+CRSNP".
	Currently, only additive effects ("+") are allowed.}
  \item{data}{object of calss \code{\link{gwaa.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{bcast}{show progress every \code{bcast} percents of progress}
  \item{simulate}{if simulated P-values should be generated}
  \item{trait.type}{Character string defining type of trait, 
	with values of "gaussian", "binomial", "poisson", "ordinal" (see help for 
	\code{\link[haplo.stats]{haplo.score.slide}} for details). If not specified, 
	the routine picks up "gaussian" or "binomial" (two levels of trait).}
  \item{\dots}{other arguments to be passed to \code{\link[haplo.stats]{haplo.score.slide}}}
}
\details{
List element P2df is set equal to P1df, as only allelic results are returned. This has 
nothing to do with actual degrees of freedom of the test.
}
\value{
	Object of class \code{\link{scan.gwaa.2D-class}}
}
\references{
For haplo.stats (scan.haplo, scan.haplo.2D), please cite:

Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland GA. (2002) Score tests for association between traits and haplotypes when linkage phase is ambiguous. Am J Hum Genet, 70: 425-434.
}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{scan.gwaa.2D-class}},
\code{\link{scan.haplo}},
\code{\link{scan.glm.2D}},
\code{\link[haplo.stats]{haplo.score.slide}}
}
\examples{
if (require(haplo.stats)) {
\dontrun{
require(GenABEL.data)
data(srdta)
c <- scan.haplo.2D("bt~sex+age+CRSNP",data=srdta,snps=(717:733),
	ids=(srdta@phdata$age<40))
plot(c)
}
}
}
\keyword{htest}% at least one, from doc/KEYWORDS
