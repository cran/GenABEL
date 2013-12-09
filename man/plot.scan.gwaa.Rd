\name{plot.scan.gwaa}
\alias{plot.scan.gwaa}
\title{function to plot GWAA results}
\description{
Plots results of GWA analysis
}
\usage{
	\method{plot}{scan.gwaa}(x, y, ..., df=1, ystart=0, col=c("blue","green"), 
		sort=TRUE, ylim, delta = 1 , main = getcall(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of type \code{\link{scan.gwaa-class}}, as returned by \code{\link{scan.glm}}, \code{\link{qtscore}},
		\code{\link{ccfast}}, \code{\link{emp.ccfast}}, \code{\link{emp.qtscore}}, or \code{\link{scan.haplo}}}
  \item{y}{this argument is not used}
  \item{\dots}{additional arguments to be passed to plot}
  \item{df}{Plot results of 1 or 2-df test (1, 2). Could be also "Pc1df" (for GC corrected P-values) 
		and "Pc2df" (for robust genomic control of the 2 df test)}
  \item{ystart}{truncate lower value of Y at this point (can help avoiding plotting too may points)}
  \item{ylim}{ylim, same as in the standard plot function}
  \item{col}{which colors to use to depict consecutive chromosomes}
  \item{sort}{whether results should be plotted after sorting by chromosome and position}
  \item{delta}{gap width between chromosomes}
  \item{main}{title of the plot}
}
%\details{}
\value{
	No value returned.
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{scan.gwaa-class}}, 
	\code{\link{add.plot}}, 
	\code{\link{snp.subset}}, 
	\code{\link{scan.glm}}, 
	\code{\link{qtscore}}, 
	\code{\link{ccfast}},
	\code{\link{emp.qtscore}}, 
	\code{\link{emp.ccfast}},
	\code{\link{scan.haplo}}
}
\examples{
require(GenABEL.data)
data(srdta)
a <- qtscore(bt,srdta,snps=c(1:250))
plot(a)
add.plot(a,df="Pc1df",col="green")
}
\keyword{hplot}% at least one, from doc/KEYWORDS
