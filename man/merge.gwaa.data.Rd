\name{merge.gwaa.data}
\alias{merge.gwaa.data}
\title{function to merge objects of gwaa.data-class}
\description{
	function to merge two objects of gwaa.data-class
}
\usage{
	\method{merge}{gwaa.data}(x, y, ...)
}
\arguments{
  \item{x}{the first object of \code{\link{gwaa.data-class}}}
  \item{y}{the second object of \code{\link{gwaa.data-class}}}
  \item{...}{arguments passed to \code{\link{merge.snp.data}}}
}
\details{
	This function calls \code{\link{merge.snp.data}} to 
	merge gtdata slots of the incoming objects; the phenotypic 
	data contained in phdata slots are merged using standard 
	function for data frames with arguments by="id" and all=TRUE. 
	The merged object is filtered and sorted according to 
	order of ids presented in the merged snp.data object.
}
\value{
	Merged object of \code{\link{gwaa.data-class}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Maksim Struchalin, Yurii Aulchenko}
%\note{ ~~further notes~~ 
%}
\seealso{
	\code{\link{merge.snp.data}},
	\code{\link{add.phdata}}
}
\examples{
	require(GenABEL.data)
	data(srdta)
	x1 <- srdta[c(1,3,5,6),c(2,4,5,6)]
	x2 <- srdta[c(2,4,5,6),c(1,3,5,6)]
	x3 <- merge(x1,x2)
	x1
	as.character(x1)
	x2
	as.character(x2)
	x3
	as.character(x3)
	srdta[1:6,1:6]
	as.character(srdta[1:6,1:6])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
