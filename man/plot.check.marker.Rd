\name{plot.check.marker}
\alias{plot.check.marker}
\title{plots "check.marker" object}
\description{
	Plots "check.marker" object, as returned by 
	\code{\link{check.marker}}
}
\usage{
	\method{plot}{check.marker}(x, y, ...)
}
\arguments{
  \item{x}{Object of class "check.marker", as returned by \code{\link{check.marker}} 
		or \code{\link{snp.subset}}}
  \item{y}{this argument is not used}
  \item{\dots}{other arguments to be passed to plot}
}
\details{
In this plot, along the X axes, you can see colour representation of 
markers which did not pass (pass -- black) the QC. The diagonal 
shows redundant markers. If for some marker there exist markers, 
which show exactly the same (or some minimum concordance) 
genotypic distribution, such 
markers are depicted as crosses an solid line is dropped on 
the X axes from it. Other solid line connects the original 
SNP with the redundant ones (depicted as circles). From each 
redundant SNP, a dashed line is dropped on X. Normally, one 
expects that redundant markers are positioned very closely and 
redundancy appears because of linkage disequilibrium. 
}
\value{
No value returned. Explanatory note is shown on the screen. 
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{check.marker}},
	\code{\link{snp.subset}}
}
\examples{
require(GenABEL.data)
data(srdta)
mc <- check.marker(data=srdta@gtdata[,1:100],redundant="all",maf=0.01,
					minconcordance=0.9,fdr=.1,ibs.mrk=0)
mc <- check.marker(data=srdta@gtdata[,1:100],maf=0.01,fdr=.1,ibs.mrk=0)
plot(mc)
mc1 <- snp.subset(mc,snps=srdta@gtdata@snpnames[20:40])
plot(mc1)
}
\keyword{hplot}
