\name{snp.mx-class}
\docType{class}
\alias{snp.mx-class}
\alias{[,snp.mx-method}
\alias{coerce,snp.mx,numeric-method}
\alias{coerce,snp.mx,character-method}
\alias{show,snp.mx-method}
\alias{summary,snp.mx-method}

\title{Class "snp.mx"}
\description{This low-level class contains objects holding large 
arrays of single nucleotide polymorphism (SNP) genotypes}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{object used to store genotypes}
	}
}
%\section{Extends}{}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "snp.mx")}: subset
      operations. [x,y] will select people listed in x and SNPs listed in y.}
    \item{coerce}{\code{signature(from = "raw", to =
	"snp.mx")}: makes an snp.mx object out of raw data}
    \item{show}{\code{signature(object = "snp.mx")}: shows the object.
      Take care that (a) this is internal representation and (b) 
	the objects are usually very large!}
	 }
}
%\references{}
\author{Yurii Aulchenko}
\note{User is not supposed to work with this class. Use \code{\link{snp.data-class}}.}
\seealso{
	\code{\link{gwaa.data-class}},
	\code{\link{snp.data-class}}
}
%\examples{}
\keyword{classes}

