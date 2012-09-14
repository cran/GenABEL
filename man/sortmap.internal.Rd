\name{sortmap.internal}
\alias{sortmap.internal}
\title{Internal function for map-sorting}
\usage{
  sortmap.internal(chrom, map, delta = 1)
}
\arguments{
  \item{chrom}{vector of markers' chromosomes}

  \item{map}{vector of marlers' map ositions}

  \item{delta}{step to do between chroms when building
  cumulative map}
}
\value{
  list, withe elements 'ix' ('sorted' order), etc.
}
\description{
  Internal function for map-sorting, not supposed to be
  used directly by user (is open for regression testing
  reasons)
}
\author{
  Yurii Aulchenko
}

