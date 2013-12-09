\name{HWE.show}
\alias{HWE.show}
\title{show HWE tables}
\description{
This function displays HWE tables and shows Chi2 and exact HWE P-values
}
\usage{
HWE.show(data, idsubset = c(1:data@gtdata@nids), 
	snpsubset = c(1:data@gtdata@nsnps))
}
\arguments{
  \item{data}{object of class \code{\link{gwaa.data-class}} or \code{\link{snp.data-class}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
}
%\details{}
\value{
	Only screen output
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{check.marker}}
}
\examples{
require(GenABEL.data)
data(srdta)
mc <- check.marker(srdta,p.lev=0.01,ibs.mrk=0)
mc$nohwe
HWE.show(data=srdta,snps=mc$nohwe)
}
\keyword{misc}% at least one, from doc/KEYWORDS
