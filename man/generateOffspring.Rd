\name{generateOffspring}
\alias{generateOffspring}
\title{simulates offspring's genotypes}
\usage{
  generateOffspring(g1, g2, q = NULL)
}
\arguments{
  \item{g1}{vector of genotypes of parent 1, coded with 0,
  1, 2 and NA}

  \item{g2}{vector of genotypes of parent 2, coded with 0,
  1, 2 and NA}

  \item{q}{vector of the frequencies of effects alleles
  (used to simulate missings in parental genotypes). If
  NULL, c(0.25,0.5,0.25) is used.}
}
\value{
  a vector containing simulated genotypes of offspring
}
\description{
  Simulat offspring's genotypes given genotypes of the
  parents. No LD is assumed.
}
\examples{
eaf <- runif(10)
g1 <- rbinom(10,2,eaf)
g2 <- rbinom(10,2,eaf)
generateOffspring(g1,g2)
}
\author{
  Yurii Aulchenko
}

