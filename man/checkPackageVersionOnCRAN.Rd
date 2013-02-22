\name{checkPackageVersionOnCRAN}
\alias{checkPackageVersionOnCRAN}
\title{checks what is the version of package on CRAN}
\usage{
  checkPackageVersionOnCRAN(packageName,
    baseUrlCRAN = "http://cran.r-project.org/web/packages/",
    timeout = 2)
}
\arguments{
  \item{packageName}{name of the package to check}

  \item{baseUrlCRAN}{path to CRAN repository}

  \item{timeout}{web chack timeout}
}
\value{
  string containing CRAN version of the package
}
\description{
  Checks what is the version of package on CRAN. The CRAN
  page (baseUrlCRAN+packageName) is checked and parsed
  extracting the line with "Package source:
  packageName_Version.tar.gz" e.g. "Package source:
  GenABEL_1.6-9.tar.gz" and then the 'Version' is returned.
  Otherwise, NULL is returned.
}
\examples{
library(GenABEL)
packageVersion("GenABEL")
checkPackageVersionOnCRAN("GenABEL")
}
\author{
  Yurii Aulchenko
}

