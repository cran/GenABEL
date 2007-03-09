\name{descriptives.trait}
\alias{descriptives.trait}
\title{Function to generate descriptive summary tables for phenotypic data}
\description{
	Function to generate descriptive summary tables for phenotypic data
}
\usage{
	descriptives.trait(data,subset,file,by=NULL,digits = 3)
}
\arguments{
  \item{data}{an object of \code{\link{snp.data-class}} or \code{\link{gwaa.data-class}}}
  \item{subset}{Subset of people to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
  \item{file}{A string specifying the name of a file to write the tables to (default is no file otput).}
  \item{by}{a binary trait; statistics will be produced seprately for the groups and compared}
  \item{digits}{number of digits to be printed}
}
\details{
}
\value{
	A table with descriptive statistics. Ptt: t-test; Pkw: kruskal.test;
	Pex: Fisher exact test (for factors with <5 levels)
}
\references{}
\author{Yurii Aulchenko}
\note{
}
\seealso{
}
\examples{
	data(srdta)
	descriptives.trait(srdta)
	descriptives.trait(srdta,by=srdta@phdata$sex)
}
\keyword{distribution}
