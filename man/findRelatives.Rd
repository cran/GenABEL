\name{findRelatives}
\alias{findRelatives}
\title{guesses relations between individuals}
\usage{
  findRelatives(gtdata, nmeivec = c(1:2), q = NULL,
    epsilon = 0.01, quiet = FALSE, OddsVsNull = 1000,
    OddsVsNextBest = 100, twoWayPenalty = log(10),
    doTwoWay = TRUE, vsIDs = NULL, gkinCutOff = NULL,
    kinshipMatrix = NULL)
}
\arguments{
  \item{gtdata}{genotypic data, either 'gwaa.data' or
  'snp.data' class, or matrix or 'databel' matrix (see
  details for format).}

  \item{nmeivec}{vector providing the degree of
  relationship to be tested (1: parent-offspring; 2: sibs,
  grandparent-grandchild; etc.).}

  \item{q}{vector of effect allele frequencis for the
  data.}

  \item{epsilon}{genotyping error rate}

  \item{quiet}{if TRUE, screen outputs supressed}

  \item{OddsVsNull}{threshold used in relationships
  inferences (see details)}

  \item{OddsVsNextBest}{threshold used in relationships
  inferences (see details)}

  \item{twoWayPenalty}{penalty on likelihoods resulting
  from models assuming two meiotic pathways}

  \item{doTwoWay}{or not}

  \item{vsIDs}{specific IDs to be tested vs others}

  \item{gkinCutOff}{if not null, sets a threshold used to
  pre-screen pairs before guessing relations. If value < 0
  provided, procedure sets threshold automatically
  (recommended)}

  \item{kinshipMatrix}{(genomic) kinship matrix (used if
  gkinCutOff!=NULL)}
}
\value{
  A list with elements call -- details of the call; profile
  -- table detailing likelihood for all pairs tested;
  estimatedNmeioses -- nids x nids matrix containing
  maximum likelihood estimate of meiotic distance for all
  pairs of individuals guess -- same as estimatedNmeioses,
  but all estimates not passing inference criteria
  (OddsVsNull, OddsVsNextBest) are NAed; compressedGuess --
  same as above, but removing cows and cols with
  missing-only elemnts
}
\description{
  This function guesses relationships (expressed as
  estimated number meiotic connection(s)) using genomic
  data. Compared to guessing relations from genomic kinship
  matrix, this procedure offers several enhancements:
}
\details{
  (1) by use of IBD/IBS 3-state space, it allows to
  distinguish between some pairs, which have the same
  kinship (e.g. parent-offspring from brother-sister;
  uncle-nephew from grandparent-grandchild, etc.)

  (2) it reports likelihood, allowing for more rigorous
  inferences

  If 'gtdata' are provided as a matrix (or 'databel'
  matrix), genotypes should be coded as 0, 1, or 2; each
  SNP corresponds to a column and each ID is a row. 'q'
  corresponds to the frequency of 'effect' (aka 'coded')
  allele, which is also equivalent to the mean(SNP)/2.0
  provided coding is correct.

  'nmeivec' is a sequence of integers, e.g. c(1,2) will
  test for parent-offspring pairs and pairs separated by
  two meioses (sibs, grandparent-grandchild, etc.). If
  'nmeivec' does not contain '0' as its first element, it
  will be automatically added (testing for twins). Also,
  nmeivec will be updated with
  c(nmeivec,max(nmeivec)+1,100) to allow for thesting of
  testing vs. 'null' (unrelated, 100) and 'most distant'
  specified by user (max(nmeivec)).

  While one may be interested to test only a sub-set of the
  data for relationships, it is recommended to provide 'q'
  estimated using all data available (see example).

  'gkinCutOff' allows use of genomic kinship matrix
  (computed internally) to pre-screen pairs to be tested.
  Use of this option with value '-1' is recommended: in
  this case threshod is set to 0.5^(max(nmeivec)+2). If not
  NULL, only pairs passing gkinCutOff are tested with the
  likelihood procedure.

  After likelihood estimation, inference on relationship is
  made. Releationship (in terms of number of meioses) is
  'guessed' if odds of likelihoods under the meiotic
  distance providing max likelihood and under the 'null'
  (maximal meiotic distance tested + 100) is greater than
  'OddsVsNull' parameter AND odds max-lik vs. the next-best
  meiotic distance is greater than 'OddsVsNextBest'
  parameter.
}
\examples{
require(GenABEL.data)
data(ge03d2.clean)
df <- ge03d2.clean[,autosomal(ge03d2.clean)]
df <- df[,sort(sample(1:nsnps(df),1000))]
eaf <- summary(gtdata(df))$"Q.2"
### donotrun
\dontrun{
relInfo <- findRelatives(df[27:30,],q=eaf)
relInfo
# look only for 1st and 2nd degree relatives
relInfo1 <- findRelatives(df[27:30],q=eaf,gkinCutOff=-1,nmeivec=c(1,2,3))
relInfo1
relInfoVS <- findRelatives(df[27:30,],q=eaf,nmeivec=c(1:6),vsIDs=idnames(df[27:30,])[1:2])
relInfoVS
}
### end norun
}
\keyword{htest}

