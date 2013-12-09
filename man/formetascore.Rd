\name{formetascore}
\alias{formetascore}
\title{function to run GWA analysis oriented for future meta-analysis}
\description{
Function to run GWA analysis -- using all functions available 
in GenABEL -- and produce output oriented for future meta-analysis
}
\usage{
formetascore(formula, data, stat = qtscore, transform = "no", 
		build = "unknown", verbosity = 1, ...)
}
\arguments{
  \item{formula}{standard formula}
  \item{data}{object of \code{\link{gwaa.data-class}}}
  \item{stat}{which GWA analysis function to apply. Could be 
	 \code{\link{mlreg}}, 
	 \code{\link{qtscore}}, 
	 \code{\link{mmscore}}, 
	 \code{\link{grammar}}, 
	 \code{\link{egscore}}, etc.
	}
  \item{transform}{Which trait transform to apply, could be
	 \code{\link{ztransform}} or 
	 \code{\link{rntransform}}. Default value is 
	"no" -- no transformation.
	}
  \item{build}{if you need that in output, specify genomic build here (e.g. "35")}
  \item{verbosity}{how much output is produced? Possible values are 0, 1, and 2}
  \item{\dots}{further arguments, passed to the ''stat'' GWA 
	analysis function}
}
\details{
	This function should be applied to analysis of 
	quantitative traits, if meta-analusis is aimed 
	afterwards.

	A transformation is applied to the formla-defined 
	residual, and the resulting trait is analysed with 
	specified function. Results are arranged as data-frame. 
}
\value{
	Data frame, containing GWA summary. The fields include:
	(1) SNP name (2) chromosome (3) position 
	(4) number of people with available data (5) effect of the allele 
	(6) standard error of the effect (7) P-value for the test 
	(8) corrected P-value (we will use Genomic Control) (9) coding, with 
	reference allele coming first (10) strand (11) frequency of the 
	reference allele (12) Exact P-value for HWE test, etc. (depends on 
	"verbosity" parameter.
}
%\references{~put references to the literature/web site here ~ }
\author{Yurii Aulchenko}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{ztransform}},\code{\link{qtscore}}}
\examples{
	require(GenABEL.data)
	data(ge03d2c)
	x <- formetascore(bmi ~ sex+age,ge03d2c)
	x[1:10,]
	x <- formetascore(bmi ~ sex+age,ge03d2c,trans=ztransform)
	x[1:10,]
	x <- formetascore(bmi ~ sex+age,ge03d2c,trans=rntransform,verbosity=2)
	x[1:10,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
