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
  The format is:
Formal class 'gwaa.data' [package "GenABEL"] with 2 slots

  ..@ phdata:`data.frame':	2500 obs. of  7 variables:

  .. ..\$ id : chr [1:2500] "p1" "p2" "p3" "p4" ...

  .. ..\$ sex: int [1:2500] 1 1 0 1 1 0 0 1 0 0 ...

  .. ..\$ age: num [1:2500] 43.4 48.2 37.9 53.8 47.5 45 52 42.5 29.7 45.8 ...

  .. ..\$ qt1: num [1:2500] -0.58 0.8 -0.52 -1.55 0.25 0.15 -0.56 NA -2.26 -1.32 ...

  .. ..\$ qt2: num [1:2500] 4.46 6.32 3.26 888 5.7 4.65 4.64 5.77 0.71 3.26 ...

  .. ..\$ qt3: num [1:2500] 1.43 3.9 5.05 3.76 2.89 1.87 2.49 2.68 1.45 0.85 ...

  .. ..\$ bt : int [1:2500] 0 1 1 1 1 0 0 1 0 0 ...

  ..@ gtdata:Formal class 'snp.data' [package "GenABEL"] with 9 slots

  .. .. ..@ nbytes    : num 625

  .. .. ..@ nids      : int 2500

  .. .. ..@ nsnps     : num 833

  .. .. ..@ idnames   : chr [1:2500] "p1" "p2" "p3" "p4" ...

  .. .. ..@ snpnames  : chr [1:833] "rs10" "rs18" "rs29" "rs65" ...

  .. .. ..@ chromosome: Named chr [1:833] "1" "1" "1" "1" ...

  .. .. .. ..- attr(*, "names")= chr [1:833] "rs10" "rs18" "rs29" "rs65" ...

  .. .. ..@ map       : Named num [1:833]  2500  3500  5750 13500 14250 ...

  .. .. .. ..- attr(*, "names")= chr [1:833] "rs10" "rs18" "rs29" "rs65" ...

  .. .. ..@ male      : int [1:2500] 1 1 0 1 1 0 0 1 0 0 ...

  .. .. ..@ gtps      :Formal class 'snp.mx' [package "GenABEL"] with 1 slots

  .. .. .. .. ..@ .Data: raw [1:625, 1:833] 55 59 55 a5 ...
}
\details{
}
\source{
}
\references{
}
\examples{
#main example: use this to see full functionality
# demo(srdta)

# load and work with srdta
data(srdta)
mc <- check.marker(srdta@gtdata)
plot(mc)
check.trait(names(srdta@phdata),srdta)
}
\keyword{datasets}
