\name{snp.data}
\alias{snp.data}
\title{creates an snp.data object}
\description{
	Creates object of class \code{\link{snp.data-class}}
}
\usage{
snp.data(nids, rawdata, idnames = as.character(c(1:nids)), 
	snpnames = as.character(c(1:(length(rawdata)/ceiling(nids/4)))), 
	chromosome = as.factor(c(1:(length(rawdata)/ceiling(nids/4)))), 
	map = as.double(c(1:(length(rawdata)/ceiling(nids/4)))), 
	male = rep(0, nids))
}
\arguments{
  \item{nids}{number of people}
  \item{idnames}{list of IDs}
  \item{male}{male indicator for IDs}
  \item{snpnames}{list of SNP names}
  \item{chromosome}{list of chromosomes SNPs belong to}
  \item{map}{map position of SNPs}
  \item{rawdata}{genotypes presented in raw data format}
}
%\details{}
\value{
  Object of class \code{\link{snp.data-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{ 
\code{\link{snp.data-class}}
}
%\examples{}
\keyword{classes}% at least one, from doc/KEYWORDS
