\name{refresh.gwaa.data}
\alias{refresh.gwaa.data}
\title{Updates an object from old to new GenABEL format}
\description{
	Attempts to update an object of gwaa.data-class 
	from old to new format
}
\usage{
refresh.gwaa.data(data,force=FALSE)
}
\arguments{
  \item{data}{An object of \code{\link{gwaa.data-class}} in pre-1.2-6 (data version 0) 
		format. 
	}
  \item{force}{When TRUE, the refreshing is forced, with any data in @strand and @coding 
		replaced by default data (0/1 coding, u-strand)
	}
}
\details{
	Takes old-style gwaa.data object and sets @coding and 
	@strand attributes to SNPs. All coding is set to 1/2 and 
	strand is set to "u" (unknown).
}
\value{
  Object of \code{\link{gwaa.data-class}} in new (GenABEL v > 1.2-6, 
  raw data format version 0.1) format.
}
%\references{
%}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{load.gwaa.data}}
}
%\examples{
%}
\keyword{manip}% at least one, from doc/KEYWORDS
