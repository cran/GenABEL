\name{getLogLikelihoodGivenRelation}
\alias{getLogLikelihoodGivenRelation}
\title{computes logLik of two blurGenotypes}
\usage{
  getLogLikelihoodGivenRelation(bGenotype1, bGenotype2,
    TransitionMatrix, q)
}
\arguments{
  \item{bGenotype1}{blurred genotype of person 1}

  \item{bGenotype2}{blurred genotype of person 2}

  \item{TransitionMatrix}{transition probability matrix}

  \item{q}{vector of effect allele frequencies}
}
\description{
  Compute logLik of genotypes of person 1 given genotypes
  of person 2 and assumed relation between the two persons
  (expressed with transition probability matrix; as
  returned with 'makeTransitionMatrix').
}
\examples{
require(GenABEL.data)
data(srdta)
# select 10 first SNPs
df <- srdta[,1:10]
# compute effect allele freq
EAF <- summary(gtdata(df))$"Q.2"
# get genotypes of first 2 people
g1 <- as.numeric(df[1:2,])
# blur all genotypes of person 1; use HWE to infer missing
bg1 <- blurGenotype(g1[1,],q=EAF)
# blur all genotypes of person 2; use HWE to infer missing
bg2 <- blurGenotype(g1[2,],q=EAF)
# generate sib-sib transision matrices
trss <- makeTransitionMatrix(EAF,nmei=c(2,2))
getLogLikelihoodGivenRelation(bg1,bg2,trss,EAF)
}
\author{
  Yurii Aulchenko
}

