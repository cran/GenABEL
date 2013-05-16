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
	YOUR\_R\_LIB\_LOCATION/exdata/srphenos.dat (phenotypes), 
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
# truncate the data to make the example faster
srdta <- srdta[seq(from=1,to=nids(srdta),by=2),seq(from=1,to=nsnps(srdta),by=2)]
mc <- check.marker(data=gtdata(srdta)[,1:100],redundant="all",
					maf=0.01,minconcordance=0.9,fdr=.1,ibs.mrk=0)
plot(mc)
check.trait(names(phdata(srdta)),srdta)
}
\keyword{datasets}
