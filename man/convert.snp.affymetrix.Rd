\name{convert.snp.affymetrix}
\alias{convert.snp.affymetrix}
\title{function to convert genotypic data from Affymetrix to internal format}
\description{
Converts genotypic data from Affymetrix format to internal genotypic data formated file
}
\usage{
convert.snp.affymetrix(dir, map, outfile, skipaffym)
}
\arguments{
  \item{dir}{
	Directory which affymetrix files storages.
	}
  \item{map}{
	File name with map (annotation) information.
	}
  \item{outfile}{
	Output data file.
	}
  \item{skipaffym}{
	Number of lines to skip in the Affymetrix file.
	}

}
\details{
	Affymetrix file has following format:

	
	some information...

	some information...

	some information...

	SNPID Call  Confidence   others column ...

	AFFX-7317060  AB  0.01709367 ...

	AFFX-7317061  BB  0.01683776 ...

	AFFX-7317067  AB  0.01704767 ...

	AFFX-7317077  AB  0.01817814 ...

	AFFX-7317078  AA  0.0006741961 ...

	AFFX-7317079  AA  0.004776776 ...

	AFFX-7317063  AB  0.006349149 ...

	AFFX-7317064  AB  0.04771883 ...

	AFFX-7317067  AA  0.04387166 ...


	The first several lines do not contain genotype information and have to be skiped. Skiped numbers of lines can be setted.
	by setting skipaffym input parameter. For above examle it has to be skipaffym=3.
	
	Every row corresponds to a SNP. The first column is snp name, the second - genotype.
	The second column can contain letters (AA, AB, BB) or figures (1, 2, 3). Another values consider as unmeasured.

	All affymetrix files must have same SNP amount and same SNP order.

	The first two lines in the map file will be skiped.

	If SNP does not exist in map (annotation) file this SNP will be skiped.
	
	Output will be writen into file pointed in outfile.

}
\value{
	Does not return any value, but writes file with GenABEL raw data
}
%\references{}
\author{Maksim Struchalin}
\note{
	The function does not check if "outfile" already exists, thus
	it is always over-written
}
\seealso{
	\code{\link{load.gwaa.data}},
	\code{\link{convert.snp.text}},
	\code{\link{convert.snp.mach}},
	\code{\link{convert.snp.tped}}
	\code{\link{convert.snp.illumina}}
}
\examples{
\dontrun{
	convert.snp.affymetrix(dir="where_is_our_aff_files", map="map_file", 
			outfile="output.raw", skipaffym=3)
}
}
\keyword{IO}
