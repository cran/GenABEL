\name{ibs}
\alias{ibs}
\title{Computes (average) Idenity-by-State for a set of people and markers}
\usage{
  ibs(data, snpsubset, idsubset = NULL,
    cross.idsubset = NULL, weight = "no", snpfreq = NULL)
}
\arguments{
  \item{data}{object of \code{\link{snp.data-class}}}

  \item{snpsubset}{Index, character or logical vector with
  subset of SNPs to run analysis on. If missing, all SNPs
  from \code{data} are used for analysis.}

  \item{idsubset}{IDs of people to be analysed. If missing,
  all people from \code{data} are used for analysis.}

  \item{cross.idsubset}{Parameter allowing parallel
  implementation. Not to be used normally. If supplied
  together with idsubset, the ibs/kinship for all pairs
  between idsubset and cross.idsubset computed.}

  \item{weight}{"no" for direct IBS computations, "freq" to
  weight by allelic frequency asuming HWE and "eVar" for
  empirical variance to be used}

  \item{snpfreq}{when option weight="freq" used, you can
  provide fixed allele frequencies}
}
\value{
  A (Npeople X Npeople) matrix giving average IBS (kinship)
  values between a pair below the diagonal and number of
  SNP genotype measured for both members of the pair above
  the diagonal.

  On the diagonal, homozygosity 0.5*(1+inbreeding) is
  provided with option 'freq'; with option 'eVar' the
  diagonal is set to 0.5; the diagonal is set to
  homozygosity with option 'no'.

  attr(computedobject,"Var") returns variance (replacing
  the diagonal when the object is used by
  \code{\link{egscore}}
}
\description{
  Given a set of SNPs, computes a matrix of average IBS for
  a group of individuals. This function facilitates quality
  control of genomic data. E.g. people with exteremly high
  (close to 1) IBS may indicate duplicated samples (or
  twins), simply high values of IBS may indicate relatives.
}
\details{
  When weight "freq" is used, IBS for a pair of people i
  and j is computed as

  \deqn{ f_{i,j} = \frac{1}{N} \Sigma_k \frac{(x_{i,k} -
  p_k) * (x_{j,k} - p_k)}{(p_k * (1 - p_k))} }

  where k changes from 1 to N = number of SNPs GW,
  \eqn{x_{i,k}} is a genotype of ith person at the kth SNP,
  coded as 0, 1/2, 1 and \eqn{p_k} is the frequency of the
  "+" allele. This apparently provides an unbiased estimate
  of the kinship coefficient.

  With "eVar" option above formula changes by using ( 2 *
  empirical variance of the genotype ) in the denominator.
  The empirical variance is computed according to the
  formula

  \deqn{ Var(g_k) = \frac{1}{M} \Sigma_i g_{ik}^2 -
  E[g_k]^2 }

  where M is the number of people

  Only with "freq" option monomorphic SNPs are regarded as
  non-informative.

  ibs() operation may be very lengthy for a large number of
  people.
}
\examples{
require(GenABEL.data)
data(ge03d2c)
set.seed(7)
# compute IBS based on a random sample of 1000 autosomal marker
selectedSnps <- sample(autosomal(ge03d2c),1000,replace=FALSE)
a <- ibs(ge03d2c,snps=selectedSnps)
a[1:5,1:5]
mds <- cmdscale(as.dist(1-a))
plot(mds)
# identify smaller cluster of outliers
km <- kmeans(mds,centers=2,nstart=1000)
cl1 <- names(which(km$cluster==1))
cl2 <- names(which(km$cluster==2))
if (length(cl1) > length(cl2)) cl1 <- cl2;
cl1
# PAINT THE OUTLIERS IN RED
points(mds[cl1,],pch=19,col="red")
# compute genomic kinship matrix to be used with e.g. polygenic, mmscore, etc
a <- ibs(ge03d2c,snps=selectedSnps,weight="freq")
a[1:5,1:5]
# now replace diagonal with EIGENSTRAT-type of diaganal to be used for egscore
diag(a) <- hom(ge03d2c[,autosomal(ge03d2c)])$Var
a[1:5,1:5]
##############################
# compare 'freq' with 'eVar'
##############################
ibsFreq <- ibs(ge03d2c,snps=selectedSnps, weight="freq")
ibsEvar <- ibs(ge03d2c,snps=selectedSnps, weight="eVar")
mdsEvar <- cmdscale( as.dist( 0.5 - ibsEvar ) )
plot(mdsEvar)
outliers <- (mdsEvar[,1]>0.1)
ibsFreq[upper.tri(ibsFreq,diag=TRUE)] <- NA
ibsEvar[upper.tri(ibsEvar,diag=TRUE)] <- NA
plot(ibsEvar,ibsFreq)
points(ibsEvar[outliers,outliers],ibsFreq[outliers,outliers],col="red")
}
\author{
  Yurii Aulchenko
}
\seealso{
  \code{\link{check.marker}},
  \code{\link{summary.snp.data}},
  \code{\link{snp.data-class}}
}
\keyword{'}
\keyword{htest}

