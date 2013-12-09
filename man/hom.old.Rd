\name{hom.old}
\alias{hom.old}
\title{function to compute average homozygosity within a person}
\description{
This function computes average homozygosity (inbreeding) for a set of 
people, across multiple markers. Can be used for Quality Control
(e.g. contamination checks)
}
\usage{
	hom.old(data, snpsubset, idsubset, weight="no")
}
\arguments{
  \item{data}{Object of \link{gwaa.data-class} or \link{snp.data-class}}
  \item{snpsubset}{Subset of SNPs to be used}
  \item{idsubset}{People for whom average homozygosity is to be computed}
  \item{weight}{When "no", homozygosity is computed as a proportion of 
		homozygous genotypes. When "freq", an estimate of 
		inbreeding coefficint is computed (see details).
		}
}
\details{
	With "freq" option, for person \eqn{i} inbreeding is estimated with 

	\deqn{
	f_i = \frac{O_i - E_i)}{(L_i - E_i)}
	}

	where \eqn{O_i} is observed homozygosity, \eqn{L_i} is the number of SNPs 
	measured in individual \eqn{i} and 

	\deqn{
	E_i = \Sigma_{j=1}^{L_i} (1 - 2 p_j (1 - p_j) \frac{T_{Aj}}{T_{Aj}-1})
	}

	where \eqn{T_{Aj}} is the numer of measured genotypes at locus \eqn{j}.

	Only polymorphic loci with number of measured henotypes >1 are used
	with this option.

	This measure is the same as used by PLINK (see reference).
	
	You should use as many people and markers as possible when estimating 
	inbreeding from marker data.
}
\value{
	With option weight="no": A matrix with rows corresponding to the ID names and colums
	showing the number of genotypes measured (NoMeasured) and 
	homozygosity (Hom).

	With option weight="freq": the same as above + expected homozygosity (E(Hom)) and
	the estimate of inbreeding, F. 
}
\references{
	Purcell S. et al, (2007) PLINK: a toolset for whole genome association and population-based 
	linkage analyses. Am. J. Hum. Genet.
}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{ibs}},
\code{\link{gwaa.data-class}},
\code{\link{snp.data-class}}
}
\examples{
require(GenABEL.data)
data(ge03d2)
h <- hom(ge03d2[,c(1:100)])
homsem <- h[,"Hom"]*(1-h[,"Hom"])/h[,"NoMeasured"]
plot(h[,"Hom"],homsem)
# wrong analysis: one should use all people (for right frequency) 
# and markers (for right F) available!
h <- hom(ge03d2[,c(1:10)])
h
}
\keyword{htest}% at least one, from doc/KEYWORDS
