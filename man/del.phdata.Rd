\name{del.phdata}
\alias{del.phdata}
\title{delete phenotypes from phdata}
\usage{
  del.phdata(data, what, all = FALSE)
}
\arguments{
  \item{data}{an object of \code{\link{gwaa.data-class}}}

  \item{what}{which phenotypes (variables) to delete,
  expressed as (vector of) names (character) or integer
  (column of phdata data frame)}

  \item{all}{if 'all'=TRUE and 'what' is misisng, all
  phenotypes are deleted, and only the 'id' and 'sex' are
  kept}
}
\description{
  This function is used to delete certain phenotypes from
  phenotypic part (phdata) of an object of
  \code{\link{gwaa.data-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
phdata(srdta)[1:5,]
srdta <- del.phdata(srdta,"qt1")
phdata(srdta)[1:5,]
srdta <- del.phdata(srdta,all=TRUE)
phdata(srdta)[1:5,]
}
\author{
  Yurii Aulchenko
}

