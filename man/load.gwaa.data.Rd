\name{load.gwaa.data}
\alias{load.gwaa.data}
\title{function to load GWAA data}
\description{
Load data (genotypes and phenotypes) from files to gwaa.data object
}
\usage{
load.gwaa.data(phenofile = "pheno.dat", genofile = "geno.raw", 
		force = TRUE, makemap = FALSE, sort = TRUE, id = "id")
}
\arguments{
  \item{phenofile}{data table with phenotypes}
  \item{genofile}{internally formatted genotypic data file (see \code{\link{convert.snp.text}} to 
			convert data)}
  \item{force}{Force loading the data if heterozygous X-chromosome genotypes are found in male}
  \item{makemap}{Make a consequtive map in case if map is provided chromosome-specifically}
  \item{sort}{Should SNPs be sorted in ascending order according to chromosome and position?}
  \item{id}{name of the column containing personal identification code in the phenofile}
}
\details{
The genofile must be the one resulting from \code{\link{convert.snp.text}},  
\code{\link{convert.snp.ped}}, \code{\link{convert.snp.tped}}, or 
\code{\link{convert.snp.illumina}} 
(see documentation for these functions for the file formats).

The phenotype file relates study subjects with their covariate and outcome values. 
In the phenotypic data file, the first line gives a description of the data contained 
in a particular column; the names should be unique, otherwise R will change them.
The first column of the phenotype file MUST contain the subjects' unique ID, named "id";
there should also be a column named "sex" and giving sex information (0 = female, 1 = male).
Other columns in the file should contain phenotypic information. Missing values should 
be coded with "NA"; binary traits should have values 0 or 1. An 
example of few first lines of a phenotype file is as follows:

id sex age bt1 qt qt1

"289982" 0 30.33 NA NA 3.93

"325286" 0 36.514 1 0.49 3.61

"357273" 1 37.811 0 1.65 5.30

"872422" 1 20.393 0 1.95 4.07

"1005389" 1 28.21 1 0.35 3.90

This file tells us that, for example, person 325286 is female (0 in second column), and she has 
"1" (usually this means a "case") value for the trait "bt1", so on. 
Person 289982 has measurements only for sex, age and qt1, while other measurements 
are missing (NA, Not Available). 

IDs are better kept in quotation (this would keep away the problem of e.g., leading zeros).
}
\value{
	Object of class gwaa.data
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{\code{\link{save.gwaa.data}}, 
\code{\link{convert.snp.text}}, 
\code{\link{convert.snp.ped}}, 
\code{\link{convert.snp.tped}},
\code{\link{convert.snp.illumina}} 
}
%\examples{}
\keyword{IO}% at least one, from doc/KEYWORDS
