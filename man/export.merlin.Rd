\name{export.merlin}
\alias{export.merlin}
\title{function to export GenABEL data in merlin format}
\description{
	Exports GenABEL data to Merlin and other pedigree formats
}
\usage{
	export.merlin(data, pedfile = "merlin.ped", datafile = "merlin.dat", 
	     mapfile = "merlin.map", format = "merlin", fixstrand = "no", 
	     extendedmap = TRUE, traits = 1, order = TRUE, stepids = 100, 
	     dpieceFun = "new")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{gwaa.data object}
  \item{pedfile}{Output pedigree data file name}
  \item{datafile}{Output data (information) file name; no output if NULL}
  \item{mapfile}{Output map file name}
  \item{format}{Output format, either "merlin" or "plink"}
  \item{fixstrand}{"no" -- the strand information and coding comes from the data; "+" -- change all coding 
	to "+" strand, "-" -- change all coding to "-" strand}
  \item{extendedmap}{if TRUE extended map (+ strand, + coding) is saved with the name "mapfile.ext",
	where "mapfile" is the parameter supplied by user}
  \item{traits}{How many fake traits to insert before first column of marker data}
  \item{order}{Should output be ordered by chromosome and position}
  \item{stepids}{make this larger for faster processing and smaller for lower RAM use}
  \item{dpieceFun}{function used to dump data. Do use default.}
}
\details{
	The use is straightforward, with only the "fixstrand" option requiring some 
	explanation. Consider a SNP on "-" strand with alleles G and A. If this SNP 
	is accessed on "+" strand, the corresponding alleles would be C and T. 
	While for example Affymetrix reports SNPs on bot "+" and "-" strands, 
	HapMap reports coding on "+" strand only. To make data compatible, and/or 
	to run imputations, one will need to convert all SNP codes to "+" strand. 
	This can be achieved by running export.merlin() with fixstrand="+" parameter.
}
\value{
	No value returned; files saved on HDD
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	To load the data to GenABEL again, use \code{\link{convert.snp.ped}},
	\code{\link{load.gwaa.data}}.
}
\examples{
#
# load(srdta)
# export.merlin(srdta[1:50,1:3])
#
}
\keyword{IO}% at least one, from doc/KEYWORDS
