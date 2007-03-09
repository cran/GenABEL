\name{GenABEL-package}
\alias{GenABEL-package}
\alias{GenABEL}
\docType{package}
\title{
GenABEL: an R package for Genome Wide Association Analysis
}
\description{

Genome-wide association (GWA) analysis is a tool of choice for identification of genes for complex traits. Effective storage, handling and analysis of GWA data represent a challenge to modern computa-tional genetics. GWA studies generate large amount of data: hun-dreds of thousands of single nucleotide polymorphisms (SNPs) are genotyped in hundreds or thousands of patients and controls. Data on each SNP undergoes several types of analysis: characterization of frequency distribution, testing of Hardy-Weinberg equilibrium, analysis of association between single SNPs and haplotypes and different traits, and so on. Because SNP genotypes in dense marker sets are correlated, significance testing in GWA analysis is pref-erably performed using computationally intensive permutation test procedures, further increasing the com-putational burden.

To make GWA analysis possible on standard desktop computers we developed GenABEL library which addresses the following objectives: 

(1) Minimisation of the amount of rapid access memory (RAM) used and the time required for data transactions. For this, we developed an effective data storage and manipulation model.

(2) Maximisation of the throughput of GWA analysis. For this, we designed optimal fast procedures for specific genetic tests. 

Imbedding GenABEL into R environment allows for easy data characterisation, exploration and presentation of the results and gives access to a wide range of standard and special statistical analysis functions available in base R and specific R packages, such as "haplo.stats", "genetics", etc.
}
\details{
\tabular{ll}{
Package: \tab GenABEL\cr
Type: \tab Package\cr
Version: \tab 1.1-8\cr
Date: \tab 2007-03-09\cr
License: \tab GNU GPL v. 2.0 or later\cr
}

To see (more or less complete) functionality of GenABEL, try running

demo(ge03d2).

Other demo of interest could be run with demo(srdta). 
Depending on your user priveleges in Windows, it may well not run. In this
case, try
demo(srdtawin).

The most important functions and classes are:

For loading the data, see
\code{\link{convert.snp.text}},
\code{\link{convert.snp.ped}},
\code{\link{load.gwaa.data}}.

Also check companion programs, affy2mega.pl and affy2gwaa.pl

For data managment, and manipulations see
\code{\link{gwaa.data-class}},
\code{\link{snp.data-class}},
\code{\link{snp.names}},
\code{\link{snp.subset}}.

For quality control, see
\code{\link{check.trait}},
\code{\link{check.marker}},
\code{\link{HWE.show}},
\code{\link{summary.snp.data}},
\code{\link{perid.summary}},
\code{\link{ibs}},
\code{\link{hom}}.

For fast analysis function, see
\code{\link{scan.gwaa-class}},
\code{\link{ccfast}},
\code{\link{qtscore}},
\code{\link{emp.ccfast}},
\code{\link{emp.qtscore}}

For specific tools facilitating analysis of the data with stratification, see
\code{\link{ibs}}. Principle components analysis (\code{\link{cmdscale}} 
is a very useful function in this context.

For functions facilitating construction of tables for your manuscript, see
\code{\link{descriptives.marker}},
\code{\link{descriptives.trait}},
\code{\link{descriptives.scan}}.

For link to WEB databases, see
\code{\link{show.ncbi}}.

For interfaces to other packages and standard R functions, 
also for 2D scans, see
\code{\link{scan.glm}},
\code{\link{scan.glm.2D}},
\code{\link{scan.haplo}},
\code{\link{scan.haplo.2D}},
\code{\link{scan.gwaa-class}},
\code{\link{scan.gwaa.2D-class}}.

For graphical facilities, see
\code{\link{plot.scan.gwaa}},
\code{\link{plot.check.marker}}.

}
\author{
Yurii Aulchenko

Maintainer: Yurii Aulchenko <i.aoultchenko@erasmusmc.nl>
}
\references{
If you use the package in your analysis, please cite the following work:

Aulchenko Y.S., Ripke S., Isaacs A., van Duijn C.M. (in press) GenABEL: an R package for genome-wide association analysis. Bioinformatics.

For exact HWE, please cite:

Wigginton G.E., Cutler D.J., Abecasis G.R. (2005) A note on exact tests of Hardy-Weinberg equilibrium. Am J Hum Genet, 76: 887-893.

For haplo.stats (scan.haplo, scan.haplo.2D), please cite:

Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland GA. (2002) Score tests for association between traits and haplotypes when linkage phase is ambiguous. Am J Hum Genet, 70: 425-434.
}
\keyword{package}
\seealso{
Packages \code{genetics},
\code{haplo.stats},
\code{qvalue},
}
\examples{
# to see complete functionality, run
#	demo(srdta)
# it will take a while!
# if demo(srdta) does not work (for Windows)m try
#	demo(srdtawin)
#
}
