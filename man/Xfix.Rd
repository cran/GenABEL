\name{Xfix}
\alias{Xfix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to set heterozygous male X-linked genotypes to NA}
\description{
Sets heterozygous male X-linked genotypes to NA
}
\usage{
	Xfix(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Object of \link{gwaa.data-class}}
}
\details{
}
\value{
	The same object of \link{gwaa.data-class}, with fixed 
	genotypes
}
\references{}
\author{Yurii Aulchenko}
\note{
}
\seealso{
\code{\link{check.marker}}
}
\examples{
data(ge03d2c)
# many errors
mc0 <- check.marker(ge03d2c)
# take only people and markers passing QC
fixed0 <- ge03d2c[mc0$idok,mc0$snpok]
# major errors fixed, still few males are heterozygous for X-chromsome markers
mc1 <- check.marker(fixed0)
# fix minor X-chromosome problems
fixed1 <- Xfix(fixed0)
# no errors
mc2 <- check.marker(fixed1)
summary(mc2)
}
\keyword{misc}% at least one, from doc/KEYWORDS
