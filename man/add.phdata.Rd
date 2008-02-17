\name{add.phdata}
\alias{add.phdata}
\title{Adds phenotypic variables to gwaa.data object}
\description{
	Adds phenotypic variables to @phdata slot 
	of an \code{\link{gwaa.data-class}} object
}
\usage{
add.phdata(data, phdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{an object of \code{\link{gwaa.data-class}}}
  \item{phdata}{data frame}
}
\details{
	The "phdata" data frame is simply merged to 
	the @phdata slot of the "data", and is sorted 
	according to the right order. 

	The "phdata" frame should contain single variable 
	named "id", preferably of "character" class. It may 
	contain "sex" variable, but that will be re-named 
	to avoid duplication with the default sex variable
	presented in @phdata.
}
\value{
	Updated object of \code{\link{gwaa.data-class}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Yurii Auclhenko}
%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
	data(srdta)
# take a small subset for this example
	srdta <- srdta[1:10,1:5]
	srdta
# generate id names
	ids <- paste("p",c(2,1,7,3,5,9,11,22,27),sep="")
# generate some random trait values
	newtra <- rnorm(9)
# make data frame
	x <- data.frame(id=ids,newtra=newtra)
	x
# now add this new trait to the data
	srdta1 <- add.phdata(srdta,x)
	srdta1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
