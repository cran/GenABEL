\name{scan.haplo}
\alias{scan.haplo}
\title{scan.haplo}
\description{
	Runs \code{\link[haplo.stats]{haplo.score.slide}} from the package 
	\code{haplo.stats} and represents output as 
	\code{\link{scan.gwaa-class}} data object
}
\usage{
scan.haplo(formula, data, snpsubset, idsubset, n.slide = 2, bcast = 10, 
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
  \item{n.slide}{Default = 2. Number of loci in each contiguous subset. The first subset is the 
		ordered loci numbered 1 to n.slide, the second subset is 2 through 
		n.slide+1 and so on. If the total number of loci in geno is n.loci, 
		then there are n.loci - n.slide + 1 total subsets.}
  \item{bcast}{show progress every \code{bcast} SNPs}
  \item{simulate}{if simulated P-values should be generated}
  \item{trait.type}{Character string defining type of trait, 
	with values of "gaussian", "binomial", "poisson", "ordinal" (see help for 
	\code{\link[haplo.stats]{haplo.score.slide}} for details). If not specified, 
	the routine picks up "gaussian" or "binomial" (two levels of trait).}
  \item{\dots}{other arguments to be passed to \code{\link[haplo.stats]{haplo.score.slide}}}
}
\details{

List element P2df is set equal to P1df, as only allelic results are returned. This has 
nothing to do with degrees of freedom.

}
\value{
	Object of class \code{\link{scan.gwaa-class}}
}
\references{
For haplo.stats (scan.haplo, scan.haplo.2D), please cite:

Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland GA. (2002) Score tests for association between traits and haplotypes when linkage phase is ambiguous. Am J Hum Genet, 70: 425-434.
}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{scan.gwaa-class}},
\code{\link[haplo.stats]{haplo.score.slide}}
}
\examples{
require(GenABEL.data)
data(srdta)
a <- ccfast("bt",srdta,snps=(717:733),ids=(srdta@phdata$age<40))
plot(a)
if (require(haplo.stats)) {
b <- scan.haplo("bt~sex+CRSNP",srdta,snps=(717:733),
		ids=(srdta@phdata$age<40))
c <- scan.haplo("bt~sex+CRSNP",srdta,snps=(717:733),
		ids=(srdta@phdata$age<40),n.slide=3)
add.plot(b,col="red",type="l")
add.plot(c,col="darkgreen",type="l")
}
}
\keyword{htest}% at least one, from doc/KEYWORDS
