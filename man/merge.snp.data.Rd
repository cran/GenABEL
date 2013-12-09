\name{merge.snp.data}
\alias{merge.snp.data}
\title{function to merge objects of snp.data-class}
\description{
	function to merge two objects of snp.data-class
}
\usage{
	\method{merge}{snp.data}(x, y, ... , error_amount = 1e+06, replacena = TRUE, 
	forcestranduse = FALSE, sort = TRUE, 
	intersected_snps_only = FALSE)
}
\arguments{
  \item{x}{the first object of \code{\link{snp.data-class}}}
  \item{y}{the second object of \code{\link{snp.data-class}}}
  \item{...}{additional arguments (not used or passed)}
  \item{error_amount}{if this amount of errors is exceeded, only error table is returned}
  \item{replacena}{Some genotypes may be missing in set 1, but measured in set 2.
	If replacena=TRUE, genotypes from the set 2 will appear in the merged data.
	}
  \item{forcestranduse}{if TRUE, forces use of strand infomation even if 
	coding information is sufficient for merging
	}
  \item{sort}{if TRUE, sorts the object according to chromosome and SNP 
	position
	}
  \item{intersected_snps_only}{if TRUE, then only interesected SNPs will be in output
	}
}
\details{
	By default, when a genotype for a person is measured in both 
	set "x" and set "y", and these are not equal, the value specified 
	by set "x" is returned in the merged set. In case when genotype 
	is NA in the first set, the behaviour depends on the value of 
	the "replacena" parameter -- if set to TRUE (default), these 
	are replaced with the non-NA values from set "y". 

	When "forcestranduse" is set to FALSE, strand information is 
	not used unless the coding is not sufficient for merging the 
	data (i.e. starnd information is used only to merge
	A/T and G/C polymorphisms). 

	SNP error is returned to "snp" table when SNP coding is 
	incompatible between the two sets. For such SNPs, only the 
	data provided by set "x" are used in the merged data.

	ID error is returned to "id" table when genotypes of the 
	same person at the same SNP are different between set "x" 
	and set "y". For such genotypes, the data provided by set "x" 
	are used in the merged data.
}
\value{
 A list is returned
  \item{id}{This table summarises individual genotype inconsistencies. 
	These may occur when some person is present and genotyped for the 
	same marker in both sets, but these genotypes are inconsistent. 
	The table's first column, "id", contains personal ID, the second, 
	"snpnames", contain SNP name, and third (set "x") and fourth 
	(set "y") contain the genotypes for this person at this SNP in 
	sets 1 and 2. 
	}
  \item{snp}{This table summarises coding errors. 
	These occur when for some SNPs coding in the set 1 is 
	not compatible with set 2. 
	The table's first column ("snpnames") provides SNP name, 
	and second (set "x") and third (set "y") report coding used in 
	respective sets.
	}
  \item{data}{merged object of \code{\link{snp.data-class}}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Maksim Struchalin, Yurii Aulchenko}
%\note{ ~~further notes~~ 
%}
\seealso{
	\code{\link{merge.gwaa.data}},
	\code{\link{add.phdata}}
}
\examples{
	require(GenABEL.data)
	data(srdta)
	x1 <- srdta[c(1,3,5,6),c(2,4,5,6)]@gtdata
	x2 <- srdta[c(2,4,5,6),c(1,3,5,6)]@gtdata
	x3 <- merge(x1,x2)
	as.character(x1)
	as.character(x2)
	as.character(x3$data)
	as.character(srdta[1:6,1:6])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
