\name{export.merlin}
\alias{export.merlin}
\title{function to export GenABEL data in merlin format}
\description{
}
\usage{
	export.merlin(data, pedfile = "merlin.ped", datafile = "merlin.dat", mapfile = "merlin.map", 
	format = "merlin", fixstrand = "no", extendedmap = TRUE, traits = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{gwaa.data object}
  \item{pedfile}{Output pedigree data file name}
  \item{datafile}{Output data (information) file name}
  \item{mapfile}{Output map file name}
  \item{format}{Output format: reserved for future use, currently only "merlin"}
  \item{fixstrand}{"no" -- the strand information and coding comes from the data; "+" -- change all coding 
	to "+" strand, "-" -- change all coding to "-" strand}
  \item{extendedmap}{if TRUE extended map (+ strand, + coding) is saved with the name "mapfile.ext",
	where "mapfile" is the parameter supplied by user}
  \item{traits}{How many fake traits to insert before first column of marker data}
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
	No value returned
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	To load the data to GenABEL again, use \code{\link{convert.snp.ped}},
	\code{\link{load.gwaa.data}}.
}
%\examples{}
\keyword{IO}% at least one, from doc/KEYWORDS
