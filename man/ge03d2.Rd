\name{ge03d2}
\alias{ge03d2}
\docType{data}
\title{GWA-type data on small region}
\description{
	\code{ge03d2} 
	A small data set (approximately 1,000 people and 8,000 SNPs) containing 
	data on 3 autosomes and X chromsome. Is a good set for 
	demonatration of the QC procedures (different genotyping errors
	are introduced) and GWA analysis. Run demo(ge03d2) to see a demo.
	This data set was developed for the "Advances in population-
	based studies" (Ge03) course of the Nihes.
}
\usage{data(ge03d2)}
\format{
}
\details{
}
\source{
}
\references{
}
\examples{
#main example: use this to see full functionality
# demo(ge03d2)

# load and work with ge03d2
data(ge03d2)
a <- qtscore(dm2,ge03d2)
plot(a)
}
\keyword{datasets}
