
\name{cocohet}
\alias{cocohet}
\alias{chi2_CG}
\title{Detecting rare recessive and compound heterozygote alleles in genome wide association.}
\description{
Detecting rare recessive and compound heterozygote alleles in genome wide association.	
}
\usage{
cocohet(data, trait, window, return_all_result=TRUE, manhettan_filename="manhettan_plot.jpeg", test="CHI2", min_expected_cut_off=-1)
}
\arguments{
  \item{data}{Genotype data for analysis. Object of class \code{\link{snp.data}}}
  \item{trait}{Vector with binary trait data. Object of class \code{\link{integer}} or \code{\link{numeric}}.}
  \item{window}{Number of SNPs on the "right" of a given SNP which are used in analysis with a SNP. Object of class \code{\link{integer}}}
  \item{return_all_result}{If FALSE then return only a vector where each element is a chisq obtained as a maximum chisq between a given SNP and SNPs on the right within a window. If TRUE then return also a matrix where chisq's for all tests are stored. Object of class \code{\link{logical}}}
  \item{manhettan_filename}{File name where manhettan plot will be saved after analysis. Object of class \code{\link{character}}}
  \item{test}{Name of the test to be performed. Object of class \code{\link{character}}}
  \item{min_expected_cut_off}{In case this is >=0 and test is NOT Pearson's chisq test then Pearson's chisq test (!) is performed only for SNPs which produce acontingency table where the expected number of subjects in each field is >min_expected_cut_off. Otherwise the specified test is performed. Object of class \code{\link{integer}} or \code{\link{numeric}}}
}
\details{

The function is an inplementation of the method aimed to detect a gene-phenotype association caused by recessive and compound heterozygote genotype states of multiple 
rare variants at a particular gene locus. The paper 'Detecting rare recessive and compound heterozygote alleles in genome wide association and sequencing studies with red hair color as example'; Fan Liu, Maksim V. Struchalin, Kate van Duijn, Albert Hofman, Andre Uitterlinden2 Yurii S. Aulchenko2 and Manfred Kayser1.
Submited to 'PLoS Genetics'.

The three tests are implemented: Pearson's chi-square test, Pearson's chi-square test with Yates correction, Fisher exact test. 
In case when the input parameter min_expected_cut_off is <0 the choosen in the input parameter "test" test is performed. If min_expected_cut_off >= 0
then always Pearson's chi-square test is
performed exept of the cases when expected number of subjects in a field of contingency table is <min_expected_cut_off. In this case the test choosen in
the input parameter test is performed.



}
\value{
 A list is returned.
  \item{chi2_max}{A vector where each element is a test statistic choosen as a maximum chisq among tests where a SNP and SNPs on the right within a window are involved.}
  \item{chi2_all}{Statistics of all tests done in the analysis. Each row of the matrix contains tests statistics for a SNP and all SNPs on the right of him within of a given window.
	For example: the ellement chi2_all[1,1] stands for a test }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Maksim Struchalin}
%\note{ ~~further notes~~ 
%}

\examples{

data(srdta)
chis2_nocorrection <- cocohet(data=srdta@gtdata, trait=srdta@phdata$bt, window=3, test="CHI2")		

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
