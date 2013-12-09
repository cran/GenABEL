\name{recodeChromosome}
\alias{recodeChromosome}
\title{Change chromosomal coding}
\usage{
  recodeChromosome(data, rules, quiet = FALSE)
}
\arguments{
  \item{data}{object of class for which 'chromosome' method
  is defined, e.g. 'gwaa.data', 'snp.data', 'scan.gwaa'}

  \item{rules}{list of pairs 'from=to'; the chromosomes
  coded in the original data set with 'from' will be
  recoded with 'to' value}

  \item{quiet}{if summary of recoding should not be printed
  to the screen}
}
\value{
  modified 'data' object
}
\description{
  Recoding of chromosomes according to the provided 'rules'
  for from -> to pairs. Most common use is anticipated when
  importing data from other software using only integers to
  represent chromosomes. In this situation the list of
  rules may look like this: list(24="X",25="Y",26="mt").
}
\note{
  'from' entries should be unique and not overlap with
  entries in 'to'
}
\examples{
require(GenABEL.data)
data(ge03d2)
table(chromosome(ge03d2))
# merge chromosome 3 and X, call chromosome 2 as 15
newdat <- recodeChromosome(ge03d2,rules=list("3"="X","2"=15))
table(chromosome(ge03d2),chromosome(newdat))
}
\author{
  Yurii Aulchenko
}
\keyword{manip}

