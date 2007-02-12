\name{ccfast}
\alias{ccfast}
\title{fast case-control analysis}
\description{
Fast case-control analysis by computing chi-square test 
from 2x2 or 2x3 tables
}
\usage{
ccfast(y, data, snpsubset, idsubset, quiet = FALSE)
}
\arguments{
  \item{y}{character name of vector of case-control status. Cases are denoted as 1 and controls as 0.}
  \item{data}{object of class "gwaa.data"}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{quiet}{do not print warning messages}
}
\details{
	This is yet not a valid function to analyse data including X chromosome. 
	Use \code{qtscore} for this purpose.
}
\value{
  Object of class \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{emp.ccfast}},
\code{\link{plot.scan.gwaa}},
\code{\link{scan.gwaa-class}}
}
\examples{
data(srdta)
a <- ccfast("bt",data=srdta,snps=c(1:10))
a
plot(a)
a <- ccfast("bt",data=srdta)
plot(a)
}
\keyword{htest}% at least one, from doc/KEYWORDS
