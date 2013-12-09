\name{mmscore}
\alias{mmscore}
\title{Score test for association in related people}
\description{
Score test for association between a trait and genetic polymorphism, 
in samples of related individuals
}
\usage{
mmscore(h2object,data,snpsubset,idsubset,strata,times=1,quiet=FALSE,
		bcast=10,clambda=TRUE,propPs=1.0)
}
\arguments{
  \item{h2object}{An object returned by \code{\link{polygenic}} polygenic mixed model analysis 
	routine. The sub-objects used are measuredIDs, residualY, and InvSigma. One can supply 
	mmscore with a fake h2object, containing these list elements. 
	}
  \item{data}{An object of \code{\link{gwaa.data-class}}. 
		ALWAYS PASS THE SAME OBJECT WHICH WAS USED FOR i\code{\link{polygenic}} ANALYSIS,
		NO SUB-SETTING IN IDs (USE IDSUBSET ARGUMENT FOR SUB-SETTING)!!!}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{strata}{Stratification variable. If provieded, scores are computed within strata and 
		then added up.}
  \item{times}{If more then one, the number of replicas to be used in derivation of 
		empirical genome-wide significance. NOTE: The structure of the data 
		is not exchangable, therefore do not use times > 1 
		unless you are really sure you understand what you are doing!
		}
  \item{quiet}{do not print warning messages}
  \item{bcast}{If the argument times > 1, progress is reported once in bcast replicas}
  \item{clambda}{If inflation facot Lambda is estimated as lower then one, this parameter 
		controls if the original P1df (clambda=TRUE) to be reported in Pc1df, 
		or the original 1df statistics is to be multiplied onto this "deflation" 
		factor (clambda=FALSE).
		If a numeric value is provided, it is used as a correction factor.}
  \item{propPs}{proportion of non-corrected P-values used to estimate the inflation factor Lambda,
		passed directly to the \code{\link{estlambda}}}
}
\details{
	Score test is performed using the formula

	\deqn{
	\frac{((G-E[G]) V^{-1} residualY)^2}{(G-E[G]) V^{-1} (G-E[G])}
	}

	where \eqn{G} is the vector of genotypes (coded 0, 1, 2) and \eqn{E[G]} is 
	a vector of (strata-specific) mean genotypic values; \eqn{V^{-1}} is the 
	InvSigma and \eqn{residualY} are residuals from the trait analysis 
	with \code{\link{polygenic}} procedure.

	This test is similar to that implemented by Abecasis et al. (see reference). 
}
\value{
  Object of class \code{\link{scan.gwaa-class}}; only 1 d.f. test is 
  implemented currently.
}
\references{
Chen WM, Abecasis GR. Family-based association tests for genome-wide association 
scans. Am J Hum Genet. 2007 Nov;81(5):913-26. 
}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{grammar}},
\code{\link{qtscore}},
\code{\link{egscore}},
\code{\link{plot.scan.gwaa}},
\code{\link{scan.gwaa-class}}
}
\examples{
# ge03d2 is rather bad data set to demonstrate, 
# because this is a population-based study
require(GenABEL.data)
data(ge03d2.clean)
#take half for speed
ge03d2.clean <- ge03d2.clean[1:100,]
gkin <- ibs(ge03d2.clean,w="freq")
h2ht <- polygenic(height ~ sex + age,kin=gkin,ge03d2.clean)
h2ht$est
mm <- mmscore(h2ht,data=ge03d2.clean)
# compute grammar
gr <- qtscore(h2ht$pgres,data=ge03d2.clean,clam=FALSE)
#compute GC
gc <- qtscore(height ~ sex + age,data=ge03d2.clean)
#compare
plot(mm,df="Pc1df",cex=0.5)
add.plot(gc,df="Pc1df",col="red")
add.plot(gr,df="Pc1df",col="lightgreen",cex=1.1)
# can see that mmscore and grammar are quite the same... in contrast to GC
}
\keyword{htest}% at least one, from doc/KEYWORDS
