\name{extract.annotation.impute}
\alias{extract.annotation.impute}
\title{extracts SNP annotation from IMPUTE files}
\usage{
  extract.annotation.impute(genofile, infofile,
    chromosome = NA, order_geno_snp_a0_a1 = c(2, 4:5),
    skip_geno = 0,
    order_info_snp_pos_freq1_info_qual_type = c(2:7),
    skip_info = 1, allow_duplicated_names = FALSE)
}
\arguments{
  \item{genofile}{IMPUTE genotype file name}

  \item{infofile}{IMPUTE info-file name}

  \item{chromosome}{chromosome}

  \item{order_geno_snp_a0_a1}{which columns to extract from
  geno-file, and what is the order for snp name, a0, and
  a0? (default is OK)}

  \item{skip_geno}{how many lines of geno-file are to be
  skept? (default is OK)}

  \item{order_info_snp_pos_freq1_info_qual_type}{which
  columns to extract from info-file, and what is the order
  for SNP name, position, frequency of allele 1, info
  (Rsq), and quality (average max post prob)? Dafult works
  for IMPUTE v2.0, but has to be changed for other
  versions. Always check!}

  \item{skip_info}{how many lines of info-file are to be
  skept before information starts? IMPUTE v2.0 has a header
  line, therefore skip_info=1 works fine; this may be
  different for other versions of IMPUTE}

  \item{allow_duplicated_names}{if duplicated SNP names are
  allowed (same order in geno and info- files is assumed
  then)}
}
\value{
  data frame containing annotaton
}
\description{
  This function extracts SNP annotation information from
  IMPUTE files. The major problem at the moment that
  info-file format of IMPUTE is a little bit unstable
  (reported information and column order varies between
  impute v1, v2, and beta-version). Therefore take special
  care to read specification of
  'order_info_snp_pos_freq1_info_qual_type'
}
\author{
  Yurii Aulchenko
}
\keyword{IO}
\keyword{manip}

