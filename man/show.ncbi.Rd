\name{show.ncbi}
\alias{show.ncbi}
\title{Shows the region on NCBI map}
\description{
This function calls web browser and direct it 
to NCBI MapViewer, to show the region of interest.
}
\usage{
show.ncbi(region)
}
\arguments{
  \item{region}{a vector containing regional landmarks} 
}
\details{
  	The elements of input vector could be SNP rs-names
}
%\value{
%}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
%\seealso{
%}
\examples{
	\dontrun{
	show.ncbi(c("rs7926624","rs11564708"))
	}
}
\keyword{misc}% at least one, from doc/KEYWORDS
