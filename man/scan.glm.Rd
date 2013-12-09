\name{scan.glm}
\alias{scan.glm}
\title{Scan GWA data using glm}
\description{
Scan GWA data using glm
}
\usage{
scan.glm(formula, family = gaussian(), data, snpsubset, idsubset, 
		bcast = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{character string containing formula to be used in \code{\link[stats]{glm}}. 
		You should put CRSNP argument in the formula,
		to arrange how the SNP from the list would be treated. This allows to 
		put in an interaction term.}
  \item{family}{family to be passed to \code{\link[stats]{glm}}}
  \item{snpsubset}{Index, character or logical vector with subset of SNPs to run analysis on. 
		If missing, all SNPs from \code{data} are used for analysis.}
  \item{idsubset}{Index, character or logical vector with subset of IDs to run analysis on. 
		If missing, all people from \code{data/cc} are used for analysis.}
  \item{data}{object of class "gwaa.data"}
  \item{bcast}{show progress every \code{bcast} SNPs}
}
%\details{
%}
\value{
  Object of class \code{\link{scan.gwaa-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{ccfast}},
	\code{\link{qtscore}},
	\code{\link{scan.gwaa-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
a <- scan.glm("bt~sex+age+CRSNP",family=binomial(),data=srdta,snps=(1:10),bcast=2)
#plot(a)

osnp <- "rs4934"
maposnp <- srdta@gtdata@map[osnp]
maposnp
reg <- snp.names(srdta,begin=maposnp-100000,end=maposnp+100000,chrom="1")
\dontrun{
a <- scan.glm("qt3~sex+age+CRSNP",data=srdta,snps=reg)
plot(a)
plot(a,df=1)
add.plot(a,df=2)
}

# interaction with sex
\dontrun{
a <- scan.glm("qt3~age+sex*CRSNP",data=srdta,snps=reg)
plot(a,df=1)
add.plot(a,df=2)
# you can do interaction with a selected polymorphisms in the same way
}
}
\keyword{htest}% at least one, from doc/KEYWORDS
