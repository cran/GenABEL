\name{cocohet}
\alias{chi2_CG}
\alias{cocohet}
\title{Test for compound heterozygote effects}
\usage{
  cocohet(data, trait, window, return_all_result = TRUE,
    makePlot = FALSE, test = "CHI2",
    min_expected_cut_off = -1)
}
\arguments{
  \item{data}{Genotype data for analysis. Object of class
  \code{\link{snp.data}}}

  \item{trait}{Vector with binary trait data. Object of
  class \code{\link{integer}} or \code{\link{numeric}}.}

  \item{window}{Number of SNPs on the "right" of a given
  SNP which are used in analysis with a SNP. Object of
  class \code{\link{integer}}}

  \item{return_all_result}{If FALSE then return only a
  vector where each element is a chisq obtained as a
  maximum chisq between a given SNP and SNPs on the right
  within a window. If TRUE then return also a matrix where
  chisq's for all tests are stored. Object of class
  \code{\link{logical}}}

  \item{makePlot}{whether the Manhattan-type plot should be
  produced (TRUE or FALSE)}

  \item{test}{Name of the test to be performed. Available
  tests are "CHI2", "YATES" (chi2 with Yates correction),
  and "FISHER". Object of class \code{\link{character}}}

  \item{min_expected_cut_off}{In case this is >=0 and test
  is NOT Pearson's chisq test then Pearson's chisq test (!)
  is performed only for SNPs which produce acontingency
  table where the expected number of subjects in each field
  is >min_expected_cut_off. Otherwise the specified test is
  performed. Object of class \code{\link{integer}} or
  \code{\link{numeric}}}
}
\value{
  A list is returned.
}
\description{
  Detecting rare recessive and compound heterozygote
  alleles in genome wide association.
}
\details{
  The function is an inplementation of the method aimed to
  detect a gene-phenotype association caused by recessive
  and compound heterozygote genotype states of multiple
  rare variants at a particular gene locus. This method is
  described in 'Detecting Low Frequent Loss-of-Function
  Alleles in Genome Wide Association Studies with Red Hair
  Color as Example'; Fan Liu, Maksim V. Struchalin, Kate
  van Duijn, Albert Hofman, Andre G. Uitterlinden, Yurii S.
  Aulchenko, and Manfred Kayser. PLoS ONE 6(11): e28145.
  doi:10.1371/journal.pone.0028145

  The three tests are implemented: Pearson's chi-square
  test, Pearson's chi-square test with Yates correction,
  Fisher exact test.  In case when the input parameter
  min_expected_cut_off is <0 the choosen in the input
  parameter "test" test is performed. If
  min_expected_cut_off >= 0 then always Pearson's
  chi-square test is performed exept of the cases when
  expected number of subjects in a field of contingency
  table is <min_expected_cut_off. In this case the test
  choosen in the input parameter test is performed.
}
\examples{
require(GenABEL.data)
data(srdta)
chis2_nocorrection <- cocohet(data=gtdata(srdta),
trait=phdata(srdta)$bt, window=3, test="CHI2")
}
\author{
  Maksim Struchalin
}
\references{
  Fan Liu, Maksim V. Struchalin, Kate van Duijn, Albert
  Hofman, Andre G. Uitterlinden, Yurii S. Aulchenko, and
  Manfred Kayser. Detecting Low Frequent Loss-of-Function
  Alleles in Genome Wide Association Studies with Red Hair
  Color as Example'. PLoS ONE 6(11): e28145.
  doi:10.1371/journal.pone.0028145
}
\keyword{manip}

