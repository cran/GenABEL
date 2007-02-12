\name{snp.subset}
\alias{snp.subset}
\title{function to subset objects of class scan.gwaa and check.marker}
\description{
Computing obejcts of class scan.gwaa may take long, especially when 
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
data(srdta)
# processing check.marker object
mc <- check.marker(srdta,snps=c(1:100),mincon=0.8)
plot(mc)
mc1 <- snp.subset(mc,snps=srdta@gtdata@snpnames[20:50])
plot(mc1)
# processing scan.gwaa object
a <- scan.glm("qt3~sex+age+CRSNP",data=srdta,snps=(1:30))
plot(a)
a1 <- snp.subset(a,snps=srdta@gtdata@snpnames[10:20])
plot(a1)
}
\keyword{misc}% at least one, from doc/KEYWORDS
