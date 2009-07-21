\name{srdta}
\alias{srdta}
\docType{data}
\title{GWA-type data on small region}
\description{
	\code{srdta} contains gwaa.data object with results 
	on a small region of 
	about 2.5 Mb. 833 SNPs are typed on 2500 people.
	NA rate is 95\%. Sex, age, two quantitative (qt1 and 
	qt2) and one binary (bt) traits are available for 
	analysis. Run demo(srdta) and check tut-srdta.pdf 
	to see examples of work with this data set. 
	Original data files used for this set are located at 
	YOUR\_R\_LIB\_LOCATION/exdata/srphenos.dat (pehnotypes), 
	srgenos.dat (human-readable genotypes) and srgenos.raw 
	(genotypes in internal format)}
\usage{data(srdta)}
\format{
	Standard object of class \code{\link{gwaa.data-class}}
}
%\details{
%}
%\source{
%}
%\references{
%}
\examples{

\dontrun{
demo(srdta)
}

# load and work with srdta
data(srdta)
mc <- check.marker(data=srdta@gtdata[,1:100],redundant="all",maf=0.01,minconcordance=0.9,fdr=.1,ibs.mrk=0)
plot(mc)
check.trait(names(srdta@phdata),srdta)
}
\keyword{datasets}
