\name{descriptives.scan}
\alias{descriptives.scan}
\title{Function to describe "top" hits in GWA scan}
\description{
	Describes "top" hits in GWA scan
}
\usage{
	descriptives.scan(data,file,top=10,sortby="P1df",digits = 6)
}
\arguments{
  \item{data}{an object of \code{\link{snp.data-class}} or \code{\link{gwaa.data-class}}}
  \item{file}{A string specifying the name of a file to write the tables to (default is no file output).}
  \item{top}{How many "top" hits to describe}
  \item{sortby}{How to pick up "top" hits ("P1df","P2df","Pgw1df","Pgw2df")}
  \item{digits}{number of digits to be printed}
}
\details{
}
\value{
	A descritive table
}
\references{}
\author{Yurii Aulchenko}
\note{
}
\seealso{
}
\examples{
	data(srdta)
	x <- qtscore(qt2~CRSNP,srdta)
	descriptives.scan(x)
}
\keyword{distribution}
