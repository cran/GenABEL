\name{reconstructNPs}
\alias{reconstructNPs}
\title{reconstruct nuclear families}
\usage{
  reconstructNPs(relationshipGuessMatrix, sex)
}
\arguments{
  \item{relationshipGuessMatrix}{meiotic relationship
  matrix, as estimated by \code{\link{findRelatives}}}

  \item{sex}{Sex, coded with 1 for males and 0 for females}
}
\value{
  A matrix containing reconstructed pedigree(s) coded in
  linkage-like format. If "fid" is zero, this means that a
  pedigree could not be assigned.
}
\description{
  Function for reconstruction of nulecar families or
  extended pedigrees based on results of
  \code{\link{findRelatives}} output (estimatd meiotic
  distance matrix). Reconstruction is based on the fact
  that parent-offspring pairs have meiotic distance of '1',
  and sibs have a distance '2+2'. If both parents and the
  offspring are genotyped, expected distace between
  offspring and both parents is '1', and distance between
  two parents is >2 (coded as 'NA').
}
\examples{
nloci <- 100
q <- runif(nloci,min=0.05,max=0.95)
# g7---g8
#   _|_
#  |   |
# g9  g10---g11
#       __|_
#      |   /\\
#    g12 g13 g14
#
nids <- 8
sex <- c(1,0,0,1,0,0,0,0)
names(sex) <- paste("g",c(7:14),sep="")
gt <- matrix(ncol=nloci,nrow=nids)
rownames(gt) <- paste("g",c(7:14),sep="")
gt["g7",] <- rbinom(nloci,2,q)
gt["g8",] <- rbinom(nloci,2,q)
gt["g11",] <- rbinom(nloci,2,q)
gt["g9",] <- generateOffspring(gt["g7",],gt["g8",],q=q)
gt["g10",] <- generateOffspring(gt["g7",],gt["g8",],q=q)
gt["g12",] <- generateOffspring(gt["g10",],gt["g11",],q=q)
gt["g13",] <- generateOffspring(gt["g10",],gt["g11",],q=q)
gt["g14",] <- gt["g13",]
aa<-findRelatives(gt,q=q,nmei=c(1:2))
aa$guess
aaPed <- reconstructNPs(aa$guess,sex)
aaPed
}
\author{
  Yurii Aulchenko
}

