\name{snp.subset}
\alias{snp.subset}
\title{function to subset objects of class scan.gwaa and check.marker}
\description{
Computing objects of class scan.gwaa may take long, especially when 
haplotypic analysis is performed. Therefore this function helps substracting 
results on some region (indicated by list of SNPs)
}
\usage{
snp.subset(data, snpsubset)
}
\arguments{
  \item{data}{object of class \code{\link{scan.gwaa-class}} or \code{\link{check.marker-class}}}
  \item{snpsubset}{character vector of snps to select}
}
%\details{}
\value{
	Object of class \code{\link{scan.gwaa-class}} or \code{\link{check.marker-class}} 
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
\code{\link{scan.gwaa-class}},
\code{\link{check.marker-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
# processing check.marker object
#mc <- check.marker(data=srdta@gtdata[,1:100],redundant="all",maf=0.01,
# minconcordance=0.9,fdr=.1,ibs.mrk=0)
mc <- check.marker(data=srdta@gtdata[,1:100],maf=0.01,fdr=.1,ibs.mrk=0)
summary(mc)
#plot(mc)
mc1 <- snp.subset(mc,snps=srdta@gtdata@snpnames[20:50])
summary(mc1)
#plot(mc1)
# processing scan.gwaa object
a <- qtscore(qt3~sex+age,data=srdta)
plot(a)
a1 <- snp.subset(a,snps=srdta@gtdata@snpnames[10:20])
plot(a1)
}
\keyword{misc}% at least one, from doc/KEYWORDS
