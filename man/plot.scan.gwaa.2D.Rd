\name{plot.scan.gwaa.2D}
\alias{plot.scan.gwaa.2D}
\title{function to plot 2D scan results}
\description{
Plots results of 2D analysis produced by
\code{\link{scan.glm.2D}} or \code{\link{scan.haplo.2D}}
}
\usage{
	\method{plot}{scan.gwaa.2D}(x, y, ..., df=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of type \code{\link{scan.gwaa.2D-class}}, as returned by \code{\link{scan.glm.2D}} or \code{\link{scan.haplo.2D}}}
  \item{y}{this argument is not used}
  \item{\dots}{additional arguments to be passed to plot}
  \item{df}{Whether 1, 2, or "all" d.f.s should be plotted. Note that for
		\code{\link{scan.haplo.2D}} 1 and 2 d.f. list the same values.
		}
}
\details{
Now plots only "allelelic" results. This is fine for \code{\link{scan.haplo.2D}}as only alleic tests are produced; however, \code{\link{scan.glm.2D}} also 
produces "genotypic" results.
}
\value{
	No value returned.
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{scan.gwaa.2D-class}}, 
	\code{\link{scan.glm.2D}}, 
	\code{\link{scan.haplo.2D}}
}
\examples{
require(GenABEL.data)
data(srdta)
a <- scan.glm.2D("qt3~CRSNP",data=srdta,snps=c(1:10))
# "allelic" results
plot(a)
# to plot "genotypic" results:
filled.contour(x=a$map,y=a$map,z=-log10(a$P2df))
}
\keyword{hplot}% at least one, from doc/KEYWORDS
