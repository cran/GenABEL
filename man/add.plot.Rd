\name{add.plot}
\alias{add.plot}
\title{function to plot additional GWAA results}
\description{
Add plot of results of GWA analysis
}
\usage{
add.plot(x, ..., df = 1,  col=c("lightgreen","lightblue"), sort=TRUE, delta = 1)
}
\arguments{
  \item{x}{object of type scan.gwaa-class, as returned by \code{\link{scan.glm}}, \code{\link{qtscore}},
		\code{\link{ccfast}}, \code{\link{emp.ccfast}}, \code{\link{emp.qtscore}}, or \code{\link{scan.haplo}}; 
	   or of type scan.gwaa.2D-class, as returned by \code{\link{scan.haplo.2D}} or \code{\link{scan.glm.2D}}.
	}
  \item{\dots}{additional arguments to be passed to plot}
  \item{df}{P-value at which df to add (1, 2 or "Pc1df")}
  \item{col}{which colors to use to depict consecutive chromosomes}
  \item{sort}{whether results should be plotted after sorting by chromosome and position}
  \item{delta}{gap width between chromosomes}
}
%\details{}
\value{
	No value returned.
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{plot}}, 
	\code{\link{snp.subset}}, 
	\code{\link{scan.glm}}, 
	\code{\link{qtscore}}, 
	\code{\link{ccfast}},
	\code{\link{emp.qtscore}}, 
	\code{\link{emp.ccfast}},
	\code{\link{scan.haplo}},
	\code{\link{scan.haplo.2D}},
	\code{\link{scan.glm.2D}}
}
\examples{
require(GenABEL.data)
data(srdta)
a <- ccfast("bt",srdta,snps=c(1:100))
plot(a)
a1 <- qtscore(bt,srdta,snps=c(1:100))
add.plot(a1,col="red",type="l")
}
\keyword{aplot}% at least one, from doc/KEYWORDS
