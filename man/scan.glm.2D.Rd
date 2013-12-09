\name{scan.glm.2D}
\alias{scan.glm.2D}
\title{Scans regional data allowing for gene-gene interaction using glm}
\description{
Scans regional data allowing for gene-gene interaction using glm
}
\usage{
scan.glm.2D(formula, family = gaussian(), data, snpsubset, idsubset, 
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
\details{
For each pair of SNPs, say snp1 and snp2, scan.glm.2D estimates 5 models. Let us denote snp1 
when it is coded as allele dose (0,1, 2) and thus results in additive model as snp1dose and 
when it is coded as 'factor' (genotypic model) as snp1factor

m00: y ~ mu [1 regression coefficient to estimate]

m10: y ~ mu + snp1dose + snp2dose [3 coefficients]

m11: y ~ mu + snp1dose + snp2dose + snp1dose * snp2dose [4 coefficients]

m20: y ~ mu + snp1factor + snp2factor [5 coefficients]

m21: y ~ mu + snp1factor + snp2factor + snp1factor * snp2factor [9 coefficients]

In the output, "P1df" refers to the test of m00 vs m10 (this is actually 2 df test); 
"P2df" refers to the test of m00 vs m20 (4 df test); "Pint1df" refers to the test of m10 
vs m11 (1 df test);  "Pint2df" refers to the test of m20 vs m21 (4 df test). The output 
is in matrix format as these P-values are generated for each pair of SNPs in turn. 
}

\value{
  Object of class \code{\link{scan.gwaa.2D-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{scan.gwaa.2D-class}},
	\code{\link{scan.haplo.2D}}
}
\examples{
\dontrun{
require(GenABEL.data)
data(srdta)
a <- scan.glm.2D("bt~sex+age+CRSNP",family=binomial(),data=srdta,snps=(1:10),bcast=2)
plot(a)
}
}
\keyword{htest}% at least one, from doc/KEYWORDS
