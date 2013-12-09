\name{Xfix}
\alias{Xfix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to set impossible genotypes as missing}
\description{
Sets impossible genotypes (e.g. heterozygous male X-linked genotypes) to missing
}
\usage{
	Xfix(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Object of \link{gwaa.data-class}}
}
\details{
	Sets to missing genotypes in the following situations:
	(1) heterozygous male X-genotypes (2) heterozygous 
	Y- and mtDNA genotypes (3) any Y-genotypes in females.
	Should only be used after \code{\link{check.marker}}, 
	which identifies systematic sex errors.
}
\value{
	The same object of \link{gwaa.data-class}, with fixed 
	genotypes
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{check.marker}}
}
\examples{
require(GenABEL.data)
data(ge03d2c)
# truncate the data to make the example faster
ge03d2c <- ge03d2c[seq(from=1,to=nids(ge03d2c),by=2),seq(from=1,to=nsnps(ge03d2c),by=2)]
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
