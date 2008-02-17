\name{ge03d2}
\alias{ge03d2}
\alias{ge03d2.clean}
\alias{ge03d2c}
\alias{ge03d2ex}
\alias{ge03d2ex.clean}
\alias{pop2.ids}
\alias{qc.nostr}
\alias{qc.str}
\docType{data}
\title{GWA-type data on few small region}
\description{
	\code{ge03d2} 
	A small data set (approximately 1,000 people and 8,000 SNPs) containing 
	data on 3 autosomes and X chromsome. Is a good set for 
	demonatration of the QC procedures (different genotyping errors
	are introduced) and GWA analysis. Run demo(ge03d2) to see a demo.
	This data set was developed for the "Advances in population-
	based studies" (Ge03) course of the Nihes.

	\code{ge03d2c} 
	A small data set (approximately 200 people and 8,000 SNPs) containing 
	data on 3 autosomes and X chromsome. This data set is 
	complementary to \code{\link{ge03d2}}.

	\code{ge03d2ex} 
	A small data set (approximately 150 people and 4,000 SNPs) containing 
	data on 3 autosomes and X chromsome. Is a good set for 
	demonatration of the QC procedures (different genotyping errors
	are introduced) and GWA analysis. 
	This data set was developed for the "Advances in population-
	based studies" (Ge03) course of the Nihes.
	See vignette "GenABEL-tutorial.pdf" for details.

	The data sets with extension ".clean" are sets after QC.
}
\usage{data(ge03d2)}
%\format{
%}
%\details{
%}
%\source{
%}
%\references{
%}
\examples{
#main example: use this to see full functionality
# demo(ge03d2)

# load and work with ge03d2
data(ge03d2)
a <- qtscore(dm2,ge03d2)
plot(a)
}
\keyword{datasets}
