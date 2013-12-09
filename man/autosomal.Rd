\name{autosomal}
\alias{autosomal}
\title{Function telling all autosomal SNPs}
\description{
	Function telling all autosomal SNPs
}
\usage{
autosomal(data)
}
\arguments{
  \item{data}{object of gwaa.data-class or snp.data-class}
}
\details{
	For every SNP, looks up the chromosome, and, when 
	it is an utosome (not X, Y, XY, or mt), reports 
	the name back
}
\value{
	Vector of SNP names
}
%\references{ ~put references to the literature/web site here ~ }
\author{Yurii Aulchenko}
%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
	require(GenABEL.data)
	data(ge03d2)
	autosomal(ge03d2)[1:10]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
