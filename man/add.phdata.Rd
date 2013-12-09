\name{add.phdata}
\alias{add.phdata}
\title{Adds phenotypic variables to gwaa.data object}
\usage{
  add.phdata(data, newphdata, name)
}
\arguments{
  \item{data}{an object of \code{\link{gwaa.data-class}}}

  \item{newphdata}{data frame or a vector with new
  phenotypic data}

  \item{name}{if 'newphdata' is a vector, the name of new
  variable should be specified in 'name'}
}
\value{
  An (updated) object of \code{\link{gwaa.data-class}}
}
\description{
  Adds phenotypic variables to \code{phdata} part of an
  \code{\link{gwaa.data-class}} object
}
\details{
  If "newphdata" is a data frame, it is simply merged to
  the \code{phdata} part of the "data", and is sorted
  according to the right order. In this case, The
  "newphdata" frame should contain single variable named
  "id", preferably of "character" class. It may contain
  "sex" variable, but that will be re-named to avoid
  duplication with the default sex variable presented in
  \code{phdata}.

  If 'newphdata' is a vector, it should be of the same
  length as the number of people in the 'data' and is
  assumed to have the same order. In this case, you also
  need to supply the name of the new phenotype via the
  'name' argument
}
\examples{
require(GenABEL.data)
data(srdta)
# take a small subset for this example
srdta <- srdta[1:10,1:5]
srdta
# add single var
rnd <- rnorm(nids(srdta))
srdta1 <- add.phdata(srdta,rnd,name="random")
srdta1
# add > 1 var
# generate id names
ids <- paste("p",c(2,1,7,3,5,9,11,22,27),sep="")
# generate some random trait values
newtra1 <- rnorm(9)
newtra2 <- rnorm(9)
# make data frame
x <- data.frame(id=ids,newtra1=newtra1,newtra2=newtra2)
x
# now add this new trait to the data
srdta1 <- add.phdata(srdta,x)
srdta1
}
\author{
  Yurii Aulchenko
}
\seealso{
  \code{\link{merge.gwaa.data}}
  \code{\link{merge.snp.data}}
}
\keyword{manip}

