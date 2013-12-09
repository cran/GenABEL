\name{descriptives.scan}
\alias{descriptives.scan}
\title{Function to describe "top" hits in GWA scan}
\description{
	Describes "top" hits in GWA scan
}
\usage{
	descriptives.scan(data,file,top=10,sortby="P1df", sep ="\t")
}
\arguments{
  \item{data}{an object of \code{\link{snp.data-class}} or \code{\link{gwaa.data-class}}}
  \item{file}{A string specifying the name of a file to write the tables to (default is no file output).}
  \item{top}{How many "top" hits to describe}
  \item{sortby}{How to pick up "top" hits ("P1df","P2df","Pgw1df","Pgw2df")}
  \item{sep}{field separator (takes effect only if file argument provided)}
}
%\details{
%}
\value{
	A descritive table
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
%\seealso{
%}
\examples{
	require(GenABEL.data)
	data(srdta)
	x <- qtscore(qt2,srdta)
	descriptives.scan(x)
}
\keyword{distribution}
