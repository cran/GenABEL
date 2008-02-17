\name{save.gwaa.data}
\alias{save.gwaa.data}
\title{function to save gwaa.data object}
\description{
	Saves GenABEL data in internal format
}
\usage{
save.gwaa.data(data, phenofile = "pheno.dat", genofile = "geno.raw", 
		human = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{gwaa.data object}
  \item{phenofile}{name of file where the phenotypes will be saved to}
  \item{genofile}{name of file where the genotypes will be saved to}
  \item{human}{if human=TRUE, saves in human-readable format (to be converted to internal format later)}
}
\details{
	When running with human=TRUE, a lot of memory (and time to complete the operation) 
	is required. Probably, this option would not work because of memory limitations 
	in a GWA scan iwth more then few hundreds of people.
	This is possible to fix; drop me a message if you need that.
}
\value{
	No value returned
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{\code{\link{load.gwaa.data}}}
%\examples{}
\keyword{IO}% at least one, from doc/KEYWORDS
