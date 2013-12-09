\name{patch_strand}
\alias{patch_strand}
\title{function to change strand}
\description{
	Changes strand in gwaa.data-class object
}
\usage{
	patch_strand(data,snpid,strand,based_on="snpnames", quiet = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{gwaa.data or snp.data object}
  \item{snpid}{vector of ids of snsp (name or position)}
  \item{strand}{vector of strands ("+","-","u")}
  \item{based_on}{either "snpnames" or "map" depending on what info is provided by snpid}
  \item{quiet}{indicates if recoding report should be directed to the screen}
}
\details{
For SNPs, as identified by 'snpid', changes strand to strand specified by 'strand'
}
\value{
	object of gwaa.data or snp.data class
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
%\seealso{
%	\code{\link{export.merlin}}.
%}
\examples{
require(GenABEL.data)
data(srdta)
as.character(srdta@gtdata@strand[1:20])
a <- patch_strand(srdta,srdta@gtdata@snpnames[1:10],rep("+",10))
as.character(a@gtdata@strand[1:20])
a <- patch_strand(srdta,srdta@gtdata@map[1:10],rep("+",10),based_on="map")
as.character(a@gtdata@strand[1:20])
}
\keyword{manip}% at least one, from doc/KEYWORDS
