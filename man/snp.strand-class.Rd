\name{snp.strand-class}
\docType{class}
\alias{snp.strand-class}
\alias{[,snp.strand,ANY,missing,missing-method}
\alias{coerce,snp.strand,character-method}
\alias{show,snp.strand-method}

\title{Class "snp.strand"}
\description{
	This class contains the strands of the 
	typed SNPs
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{nucleotide strand data}
	}
}
%\section{Extends}{}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "snp.strand", i = "ANY", j = "missing", drop = "missing")}: subset
      operations. x[i] will show strand for SNPs selected in i.}
    \item{coerce}{\code{signature(from = "snp.strand", to =
	"character")}: converts SNP strand from internal (raw) to human-readable character.}
    \item{show}{\code{signature(object = "snp.strand")}: shows the object.
      Take care that this is internal representation
		}
	 }
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{snp.coding-class}},
	\code{\link{gwaa.data-class}},
	\code{\link{snp.data-class}}
}
\examples{
	require(GenABEL.data)
	data(srdta)
	srdta@gtdata@strand[1:10]
	as.character(srdta@gtdata@strand[1:10])
}
\keyword{classes}

