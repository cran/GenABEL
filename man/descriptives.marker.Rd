\name{descriptives.marker}
\alias{descriptives.marker}
\title{Function to generate descriptive summary tables for genotypic data}
\description{
	Function to generate descriptive summary tables for genotypic data
}
\usage{
	descriptives.marker(data,snpsubset,idsubset,file,mafc,hwec,snpc,idcc,digits = 3)
}
\arguments{
  \item{data}{an object of \code{\link{snp.data-class}} or \code{\link{gwaa.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data} are used for analysis.}
  \item{file}{A string specifying the name of a file to write the tables to (default is missing).}
  \item{mafc}{vector containing desired cut-off levels for minor allele frequency}
  \item{hwec}{vector containing desired cut-off levels for exact HWE P-values}
  \item{snpc}{vector containing desired cut-off levels for SNP call rate}
  \item{idcc}{vector containing desired cut-off levels for individual SNP call rate}
  \item{digits}{number of digits to be printed}
}
%\details{
%}
\value{
	A list containing descriptive tables and statistics
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
	descriptives.marker(srdta)
}
\keyword{distribution}
