\name{hom}
\alias{hom}
\title{function to compute average homozygosity within a person}
\usage{
  hom(data, snpsubset, idsubset, snpfreq, n.snpfreq = 1000)
}
\arguments{
  \item{data}{Object of \link{gwaa.data-class} or
  \link{snp.data-class}}

  \item{snpsubset}{Subset of SNPs to be used}

  \item{idsubset}{People for whom average homozygosity is
  to be computed}

  \item{snpfreq}{when option weight="freq" used, you can
  provide fixed allele frequencies}

  \item{n.snpfreq}{when option weight="freq" used, you can
  provide a vector supplying the number of people used to
  estimate allele frequencies at the particular marker, or
  a fixed number}
}
\value{
  A matrix with rows corresponding to the ID names and
  columns showing the number of SNPs measured in this
  person (NoMeasured), the number of measured polymorphic
  SNPs (NoPoly), homozygosity (Hom), expected homozygosity
  (E(Hom)), variance, and the estimate of inbreeding, F.
}
\description{
  This function computes average homozygosity (inbreeding)
  for a set of people, across multiple markers. Can be used
  for Quality Control (e.g. contamination checks)
}
\details{
  Homozygosity is measured as proportion of homozygous
  genotypes observed in a person.

  Inbreeding for person \eqn{i} is estimated with

  \deqn{ }{ f_i = ((O_i - E_i))/((L_i - E_i)) }\deqn{ f_i =
  \frac{(O_i - E_i)}{(L_i - E_i)} }{ f_i = ((O_i -
  E_i))/((L_i - E_i)) }\deqn{ }{ f_i = ((O_i - E_i))/((L_i
  - E_i)) }

  where \eqn{O_i} is observed homozygosity, \eqn{L_i} is
  the number of SNPs measured in individual \eqn{i} and

  \deqn{ }{ E_i = Sigma_(j=1)^(L_i) (1 - 2 p_j (1 - p_j)
  (T_(Aj))/(T_(Aj)-1)) }\deqn{ E_i = \Sigma_{j=1}^{L_i} (1
  - 2 p_j (1 - p_j) \frac{T_{Aj}}{T_{Aj}-1}) }{ E_i =
  Sigma_(j=1)^(L_i) (1 - 2 p_j (1 - p_j)
  (T_(Aj))/(T_(Aj)-1)) }\deqn{ }{ E_i = Sigma_(j=1)^(L_i)
  (1 - 2 p_j (1 - p_j) (T_(Aj))/(T_(Aj)-1)) }

  where \eqn{T_{Aj}} is the number of measured genotypes at
  locus \eqn{j}; \eqn{T_{Aj}} is either estimated from data
  or provided by "n.snpfreq" parameter (vector). Allelic
  frequencies are either estimated from data or provided by
  the "snpfreq" vector.

  This measure is the same as used by PLINK (see
  reference).

  The variance (Var) is estimated as

  \deqn{ V_{i} = \frac{1}{N} \Sigma_k \frac{(x_{i,k} -
  p_k)^2}{(p_k * (1 - p_k))} }

  where k changes from 1 to N = number of SNPs,
  \eqn{x_{i,k}} is a genotype of ith person at the kth SNP,
  coded as 0, 1/2, 1 and \eqn{p_k} is the frequency of the
  "+" allele.

  Only polymorphic loci with number of measured genotypes
  >1 are used with this option.

  This variance is used as diagonal of the genomic kinship
  matrix when using EIGENSTRAT method.

  You should use as many people and markers as possible
  when estimating inbreeding/variance from marker data.
}
\examples{
require(GenABEL.data)
data(ge03d2)
h <- hom(ge03d2[,c(1:100)])
h[1:5,]
homsem <- h[,"Hom"]*(1-h[,"Hom"])/h[,"NoMeasured"]
plot(h[,"Hom"],homsem)
# wrong analysis: one should use all people (for right frequency)
# and markers (for right F) available!
h <- hom(ge03d2[,c(1:10)])
h
}
\author{
  Yurii Aulchenko, partly based on code by John Barnard
}
\references{
  Purcell S. et al, (2007) PLINK: a toolset for whole
  genome association and population-based linkage analyses.
  Am. J. Hum. Genet.
}
\seealso{
  \code{\link{ibs}}, \code{\link{gwaa.data-class}},
  \code{\link{snp.data-class}}
}
\keyword{htest}

