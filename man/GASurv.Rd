\name{GASurv}
\alias{GASurv}
\title{Makes survival data object for reg.gwaa}
\description{
Helper to \code{\link{mlreg}}: makes survival data object
}
\usage{
GASurv(fuptime, status)
}
\arguments{
  \item{fuptime}{Follow-up time}
  \item{status}{status (1=event, 0=censored)}
}
%\details{}
\value{
	Matrix with column 1 = follow-up time, and 2 = status
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{mlreg}}
}
%\examples{}
\keyword{manip}
