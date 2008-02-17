\name{grammar}
\alias{grammar}
\title{Approximate score test for association in related people}
\description{
Fast approximate score test for association between a trait and genetic polymorphism, 
in samples of related individuals. When used with argument "times=1", it is equivalent to 
running \code{\link{qtscore}} on "environmental residuals" from \code{\link{polygenic}}.
However, it does not produce correct results with permutations, because the raw 
trait values, which are not exchangeable, are permuted. Use \code{\link{qtscore}} 
on "environmental residuals" when you want to have empirical GW significance with GRAMMAR 
method.
}
\usage{
grammar(h2object,data,snpsubset,idsubset,strata,times=1,quiet=FALSE,bcast=10,clambda=FALSE,propPs=1.0)
}
\arguments{
  \item{h2object}{An object returned by \code{\link{polygenic}} polygenic mixed model analysis 
	routine. The sub-objects used are measuredIDs, residualY, h2an\$estimates (last element, 
	total variance, only), and InvSigma. 
	One can supply grammar with a fake h2object, containing these list elements. 
	}
  \item{data}{An object of \code{\link{gwaa.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{strata}{Stratification variable. If provieded, scores are computed within strata and 
		then added up.}
  \item{times}{If more then one, the number of replicas to be used in derivation of 
		empirical genome-wide significance. NOTE: do not use times > 1 
		unless you are really sure you understand what you are doing!
		}
  \item{quiet}{do not print warning messages}
  \item{bcast}{If the argument times > 1, progress is reported once in bcast replicas}
  \item{clambda}{If inflation facor Lambda is estimated as lower then one, this parameter 
		controls if the original P1df (clambda=TRUE) to be reported in Pc1df, 
		or the original 1df statistics is to be multiplied onto this "deflation" 
		factor (clambda=FALSE). With GRAMMAR, Lambda is expected ot be less than 1.
		If a numeric value is provided, it is used as a correction factor.}
  \item{propPs}{proportion of non-corrected P-values used to estimate the inflation factor Lambda,
		passed directly to the \code{\link{estlambda}}}
}
\details{
	Approximate score test is performed using the formula

	$$
	\sigma^4 \frac{((G-E[G]) V^{-1} residualY)^2}{(G-E[G]) (G-E[G])}
	$$

	where $sigma^4$ is the square of the residual 
	variance, $G$ is the vector of genotypes (coded 0, 1, 2) and $E[G]$ is 
	a vector of (strata-specific) mean genotypic values; $V^{-1}$ is the 
	InvSigma and $residualY$ are residuals from the trait analysis 
	with \code{\link{polygenic}} procedure.

	Compared to score test implemented in \code{\link{mmscore}}, grammar 
	test is faster and computation time grows only linearly with the 
	number of subjects (with \code{\link{mmscore}} this relation is 
	quadratic). While raw P1df from grammar are not quite correct, 
	the GC p-values correspond very closely to these from the \code{\link{mmscore}}.
}
\value{
  Object of class \code{\link{scan.gwaa-class}}; only 1 d.f. test is 
  implemented currently.
}
\references{
Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid association using mixed model 
and regression: a fast and simple method for genome-wide pedigree-based quantitative 
trait loci association analysis. Genetics. 2007 177(1):577-85.

Amin N, van Duijn CM, Aulchenko YS. A genomic background based method for 
association analysis in related individuals. PLoS ONE. 2007 Dec 5;2(12):e1274. 
}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{grammar}},
\code{\link{qtscore}},
\code{\link{plot.scan.gwaa}},
\code{\link{scan.gwaa-class}}
}
\examples{
# ge03d2.clean is rather bad data set to demonstrate grammar, 
# because this is a population-based study
data(ge03d2.clean)
#take half for speed
ge03d2.clean <- ge03d2.clean[1:450,]
gkin <- ibs(ge03d2.clean,w="freq")
h2ht <- polygenic(height ~ sex + age,kin=gkin,ge03d2.clean)
h2ht$est
# can use "grammar", but ...
a <- grammar(h2ht,data=ge03d2.clean)
# ... use rather qtscore (note clam=FALSE), which is a better alternative for the same thing
a <- qtscore(h2ht$pgres,data=ge03d2.clean,clam=FALSE)
# compare to GC:
b <- qtscore(height ~ sex + age,data=ge03d2.clean)
plot(b,df="Pc1df",col="red")
add.plot(a,df="Pc1df")
# relatively large difference is due to high heritability
# note that locus at chromosome 2 should indeed be there...
}
\keyword{htest}% at least one, from doc/KEYWORDS
