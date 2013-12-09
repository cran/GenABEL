\name{snp.coding-class}
\docType{class}
\alias{snp.coding-class}
\alias{[,snp.coding,ANY,missing,missing-method}
\alias{coerce,snp.coding,character-method}
\alias{show,snp.coding-method}

\title{Class "snp.coding"}
\description{
	This class contains the actual nucleotide codes for the 
	typed SNPs
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{nucleotide coding data}
	}
}
%\section{Extends}{}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "snp.coding", i = "ANY", j = "missing", drop = "missing")}: subset
      operations. x[i] will show coding for SNPs selected in i.}
    \item{coerce}{\code{signature(from = "snp.coding", to =
	"character")}: converts SNP coding from internal (raw) to human-readable character.}
    \item{show}{\code{signature(object = "snp.coding")}: shows the object.
      Take care that this is internal representation
		}
	 }
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{snp.strand-class}},
	\code{\link{gwaa.data-class}},
	\code{\link{snp.data-class}}
}
\examples{
	require(GenABEL.data)
	data(srdta)
	srdta@gtdata@coding[1:10]
	as.character(srdta@gtdata@coding[1:10])
}
\keyword{classes}

