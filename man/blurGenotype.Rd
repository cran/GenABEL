\name{blurGenotype}
\alias{blurGenotype}
\title{blur genotype calls into probabilites}
\usage{
  blurGenotype(g, q = NULL, epsilon = 0.01)
}
\arguments{
  \item{g}{vector of genotypes for a particular person (at
  locus 1, locus 2, etc., coded as 0, 1, 2 (corresponding
  to genotypes AA, AB, and BB, respectively) and NA.}

  \item{q}{(optional) vector of coded allele freqeuncies
  for locus 1, locus 2, etc.}

  \item{epsilon}{error rate}
}
\value{
  matrix with columns corresponding to SNPs and rows
  corresponding to 'g0', 'g1', 'g2'. For a particular SNP,
  a vale in cell 'gK' is the probability that true genotype
  is 'K', given thw original call and error-rate.
}
\description{
  'blurs' genotype calls into probabilities: translates
  single genotype g2, into probability distribution
  P(g1|g2), that is probability that true genotype is g1
  given g2 is the observed 'called' genotype and error rate
  is epsilon. Probability that 'true' genotype is called
  genotype is set to (1-epsilon)^2, the probability that
  true genotype differs at 1 allele is set to
  epsilon*(1-epsilon), and both allelel differ = epsilon^2.
}
\examples{
require(GenABEL.data)
data(srdta)
# select 10 first SNPs
df <- srdta[,1:10]
# compute effect allele freq
EAF <- summary(gtdata(df))$"Q.2"
EAF
# get genotypes of first 5 people
g1 <- as.numeric(df[1:5,])
g1
# blur the genotype of person 1, snp 1
blurGenotype(g1[1,1])
# blur all genotypes of person 2; assume no info for missing
blurGenotype(g1[2,])
# blur all genotypes of person 2; use HWE to infer missing
blurGenotype(g1[2,],q=EAF)
}
\author{
  Yurii Aulchenko
}

