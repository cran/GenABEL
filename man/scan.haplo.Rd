\name{scan.haplo}
\alias{scan.haplo}
\title{scan.haplo}
\description{
	Runs \code{\link[haplo.stats]{haplo.score.slide}} from the package 
	\code{haplo.stats} and represents output as 
	\code{\link{scan.gwaa-class}} data object
}
\usage{
scan.haplo(y, data, snpsubset, idsubset, n.slide = 2, bcast = 25, simulate=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{character name of the vector of trait values}
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
data(srdta)
a <- ccfast("bt",srdta,snps=(717:733),ids=(srdta@phdata$age<40))
b <- scan.haplo("bt",srdta,snps=(717:733),
		ids=(srdta@phdata$age<40),trait.type="binomial")
c <- scan.haplo("bt",srdta,snps=(717:733),
		ids=(srdta@phdata$age<40),trait.type="binomial",n.slide=3)
plot(a)
add.plot.scan.gwaa(b,col="red",type="l")
add.plot.scan.gwaa(c,col="darkgreen",type="l")
}
\keyword{htest}% at least one, from doc/KEYWORDS
