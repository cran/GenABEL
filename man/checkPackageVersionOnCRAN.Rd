\name{checkPackageVersionOnCRAN}
\alias{checkPackageVersionOnCRAN}
\title{checks what is the version of package on CRAN...}
\usage{checkPackageVersionOnCRAN(packageName,
    baseUrlCRAN="http://cran.r-project.org/web/packages/", timeout=10)}
\description{checks what is the version of package on CRAN}
\details{Checks what is the version of package on CRAN.
The CRAN page (baseUrlCRAN+packageName) is checked 
and parsed extracting the line with
"Package source:	 packageName_Version.tar.gz" 
e.g. 
"Package source:	 GenABEL_1.6-9.tar.gz"
and then the 'Version' is returned. 
Otherwise, NULL is returned.}
\value{string containing CRAN version 
of the package}
\author{Yurii Aulchenko}
\arguments{\item{packageName}{name of the package to check}
\item{baseUrlCRAN}{path to CRAN repository}
\item{timeout}{web chack timeout}}
\examples{library(GenABEL)
packageVersion("GenABEL")
checkPackageVersionOnCRAN("GenABEL")}
