\name{makeTransitionMatrix}
\alias{makeTransitionMatrix}
\title{Genotype transition probabilities matrices}
\usage{
  makeTransitionMatrix(q, nmeioses = 1000)
}
\arguments{
  \item{nmeioses}{number of meioses separating two
  individuals ((a vector) of non-negative integers). If a
  vector, it is assumed it lists all meiotic paths
  connecting the pair}

  \item{q}{(a vector of) the coded allele frequency(ies)
  (e.g. "Q.2" of GenABEL-package)}
}
\value{
  If q is scalar, a 3x3 matrix is returned, where elements
  represent conditional transition probabilities
  P(g1|g2,nmeioses); rows correspond to the genotypes of
  g1, and columns correspond to the genotypes of g2. If
  coded allele is 'B', then e.g. element [1,2] gives the
  probability P(g1='AA'|g2='AB',nmeioses=nmeioses).

  If q is a vector, series of above-described matrices are
  returned as an 'array' object. A matrix constructed for
  certain element q[i] can be accessed via result[,,i].
}
\description{
  Function to generate genotypic transition probablilites
  matrices, which represent conditional probabilities
  P(g1|g2,nmeioses), where g1 is henotype of person one
  (AA, AB or BB), g2 is genotype of person two, and
  nmeioses is the number of meioses separating these two
  individuals (0 for twins, 1 for parent-offspring, c(2,2)
  for sibs, 2 for grandparent-grandchild pairs, etc.)
}
\examples{
# transition matrix for parent-offspring, for q=0.1
makeTransitionMatrix(0.1,nmeioses=1)
# for a set of q's
makeTransitionMatrix(c(0.1,0.9),nmeioses=1)
# for sibs
makeTransitionMatrix(0.1,nmeioses=c(2,2))
# for half-sibs (or grandparent-grandchild)
makeTransitionMatrix(0.1,nmeioses=2)
# for remote relatives
makeTransitionMatrix(0.1,nmeioses=10)
# for independent
makeTransitionMatrix(0.1,nmeioses=1000)
}
\author{
  Yurii Aulchenko
}

