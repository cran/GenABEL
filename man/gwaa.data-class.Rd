\name{gwaa.data-class}
\docType{class}
\alias{gwaa.data-class}
\alias{show,gwaa.data-method}
\alias{[,gwaa.data,ANY,ANY,ANY-method}
\alias{summary,gwaa.data-method}
\alias{gtdata,gwaa.data-method}
\alias{gtdata}
\alias{phdata,gwaa.data-method}
\alias{phdata}
\alias{phdata<-,gwaa.data-method}
\alias{phdata<-}

\title{Class "gwaa.data"}
\description{This class contains objects holding all GWAA data -- 
phenotypes, genotypes and other relevant information
}
\section{Slots}{
	 \describe{
    \item{\code{phdata}:}{dataframe with phenotypic data used in GWAA}
    \item{\code{gtdata}:}{object of class \code{\link{snp.data-class}} used to store genotypic data, map, etc.}
	}
}
%\section{Extends}{
%}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "gwaa.data", i = "ANY", j = "ANY", drop = "ANY")}: subset
      operations. x[i,j] will select people listed in i and SNPs listed in j.}
    \item{show}{\code{signature(object = "gwaa.data")}: shows both parts of the object.
      Take care that the objects are usually very large!}
    \item{summary}{\code{signature(object = "gwaa.data")}: Calls 
	standard summary to describe phenotypic part and calls 
	\code{\link{summary.snp.data}} to \code{\link{snp.data-class}} }
    \item{gtdata}{\code{signature(object = "gwaa.data")}: returns genotypic data 
    (object of \code{\link{snp.data-class}})}
    \item{phdata}{\code{signature(object = "gwaa.data")}: returns phenotypic data 
    (a \code{data.frame})}
    \item{phdata<-}{\code{signature(object = "gwaa.data")}: method do modify phenotypic data 
    (a \code{data.frame})}
    }
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{snp.data-class}},
	\code{\link{load.gwaa.data}},
	\code{\link{snp.mx-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
srdta@phdata[1:10,]
srdta@gtdata[1:10,1:12]
srdta[1:10,1:12]
as.numeric(srdta@gtdata[1:12,1:10])
# very long output:
summary(srdta)
}
\keyword{classes}

