\name{check.marker}
\alias{check.marker}
\title{function to do genotypic quality control}
\description{
This function helps selecting the marker which should enter into 
GWA analysis based on call rate, minor allele frequency, 
value of the chi-square test for Hardy-Weinberg equilibrium, 
and redudndance, defined as concordance between the distributions 
of the genotypes (including missing values).
}
\usage{
check.marker(data, snpsubset, idsubset, callrate = 0.95, 
     perid.call=0.95, extr.call = 0.1, extr.perid.call = 0.1, het.fdr = 0.01, 
     ibs.threshold = 0.95, ibs.mrk = 2000, ibs.exclude="both", maf, p.level = -1, 
     fdrate = 0.2, odds = 1000, hweidsubset, redundant = "no", 
     minconcordance = 2.0, qoption = "bh95", imphetasmissing = TRUE, XXY.call=0.8,
     intermediateXF = c(0.5, 0.5)) 
}
\arguments{
  \item{data}{gwaa.data or snp.data object}
  \item{snpsubset}{a subset of SNPs to check (names, indexes, logical), default is all from \code{data}}
  \item{idsubset}{a subset of people to check (names, indexes, logical), default is all from \code{data}}
  \item{callrate}{cut-off SNP call rate}
  \item{perid.call}{cut-off individual call rate (maximum percent of missing genotypes in a person)}
  \item{extr.call}{SNPs with this low call rate are dropped prior to main analysis}
  \item{extr.perid.call}{people with this low call rate are dropped prior to main analysis}
  \item{het.fdr}{FDR rate for unacceptably high individual heterozygosity}
  \item{ibs.threshold}{threshold value for acceptable IBS}
  \item{ibs.mrk}{How many random markers should be used to estimate IBS. When ibs.mrk < 1, IBS checks are turned off.
		When "all" all markers are used.
		}
  \item{ibs.exclude}{"both", "lower" or "none" -- whether both samples with IBS>ibs.threshold 
  		should be excluded, the one with lower call rate, or no check (equivalent to use of 
  		'ibs.mrk = -1').
		}
  \item{maf}{cut-off Minor Allele Frequency. If not specified, the default value is 5 chromosomes 5/(2*nids(data))}
  \item{p.level}{cut-off p-value in check for Hardy-Weinberg Equilibrium. If negative, FDR is applied}
  \item{fdrate}{cut-off FDR level in check for Hardy-Weinberg Equilibrium}
  \item{odds}{cut-off odds to decide whether marker/person should be excluded based on sex/X-linked marker data inconsistency}
  \item{hweidsubset}{a subset of people to check (names, indexes, logical) to use for HWE check}
  \item{redundant}{if "bychrom", redundancy is checked within chromosomes; "all" -- all pairs of markers; "no" -- no redundancy checks}
  \item{minconcordance}{a parameter passed to "redundant" function. If "minconcordance" is > 1.0 
			only pairs of markers which are exactly the same, including NA pattern,
			are considered as redundant; if 0 < "minconcordance" < 1, then pairs 
			of markers with concordance > "minconcordance" are considered redundant.
			see \code{\link{redundant}} for details. Note that if "minconcordance" < 1
			the program will take much longer time to run}
  \item{qoption}{if "bh95", BH95 FDR used; if "storey", qvalue package (if installed) is used}
  \item{imphetasmissing}{If "impossible heterozygotes" (e.g. heterozygous 
			mtDNA, and male Y- and X-chromsome markers) should 
			be treated as missing genotypes in the QC procedure}
  \item{XXY.call}{What proportion of Y-chromosome markers should be called to consider that 
			Y-chromosome is present (in presence of XX)}
  \item{intermediateXF}{X-chromosomal F to be considered 'intermediate' and regarded as error; 
  			currently use of default disables this check}
}
\details{
	In this procedure, sex errors are identified initally and then 
	possible residual errors are removed iteratively. 
	At the first step, of the iterative procedure, 
	per-marker (minor allele frequency, call rate,
	exact P-value for Hardy-Weinberg equilibrium) and between-marker 
	statistics are generated and controlled for, mostly using the 
	internal call to the function \code{\link{summary.snp.data}}.
	
	At the second step of the iterative procedure, 
	per-person statistics, such call rate within 
	a person, heterozygosity and and between-person statistics 
	(identity by state across a random sample of markers) are generated, 
	using \code{\link{perid.summary}} and \code{\link{ibs}} functions. 
	Heterozygosity and IBS are estimated using only autosomal data. 
	If IBS is over ibs.threshold for a pair, one person from the 
	pair is added to the ibsfail list and excluded from the idok list. 
	At the second step, only the markers passing the first step are used. 

	The procedure is applied recursively till no further markers and 
	people are eliminated.
%	For IBS computations, only people passing personla call rate and 
%	heterozygosity checks are used. 
}
\value{
  Object of class \code{\link{check.marker-class}}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
	\code{\link{check.trait}},
	\code{\link{ibs}},
	\code{\link{summary.snp.data}},
	\code{\link{perid.summary}},
	\code{\link{plot.check.marker}},
	\code{\link{summary.check.marker}},
	\code{\link{redundant}},
	\code{\link{HWE.show}},
	\code{\link{check.marker-class}}
}
\examples{
# usual way
require(GenABEL.data)
data(ge03d2c)
# truncate the data to make the example faster
ge03d2c <- ge03d2c[seq(from=1,to=nids(ge03d2c),by=2),seq(from=1,to=nsnps(ge03d2c),by=2)]
# many errors
mc0 <- check.marker(ge03d2c)
# take only people and markers passing QC
fixed0 <- ge03d2c[mc0$idok,mc0$snpok]
# major errors fixed, still few males are heterozygous for X-chromsome markers
mc1 <- check.marker(fixed0)
# fix minor X-chromosome problems
fixed1 <- Xfix(fixed0)
# no errors
mc2 <- check.marker(fixed1)
summary(mc2)
# ready to use fixed1 for analysis

# let us look into redundancy
require(GenABEL.data)
data(srdta)
mc <- check.marker(data=srdta,ids=c(1:300),call=.92,perid.call=.92)
names(mc)
mc$nohwe
mc <- check.marker(data=srdta@gtdata[,1:100],call=0.95,perid.call=0.9,
		maf=0.02,minconcordance=0.9,fdr=0.1,redundant="all",ibs.mrk=0)
summary(mc)
HWE.show(data=srdta,snps=mc$nohwe)
plot(mc)

}
\keyword{misc}% at least one, from doc/KEYWORDS
