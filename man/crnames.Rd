\name{crnames}
\alias{crnames}
\docType{data}
\title{Return column and row names}
\description{
  Given a dimnames, returns column and row names for index 
  cells
}
\usage{
crnames(dnames,idx)
}
\arguments{
	\item{dnames}{object \code{\link{dimnames}}}
	\item{idx}{index (or logical condition on the original object)}
}
%\details{
%}
%\source{
%}
%\references{
%}
\examples{
	require(GenABEL.data)
	data(ge03d2ex)
	a <- as.numeric(ge03d2ex[1:20,1:3])
	crnames(dimnames(a),a==1)
}
\keyword{misc}
