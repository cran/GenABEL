#ifndef __export_plink_H__
#define __export_plink_H__

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif

void get_snps_many(char *a, int *Nsnps, int *Nrows, int *b);

std::string* getGenotype(std::string coding, std::string sep);

SEXP export_plink(SEXP idnames, SEXP snpdata, SEXP Nsnps, SEXP NidsTotal, SEXP Coding, SEXP from, SEXP to,
		SEXP male, SEXP traits, SEXP pedfilename, SEXP plink, SEXP append);

SEXP export_plink_tped(SEXP snpnames, SEXP chromosomes, SEXP map,
		SEXP snpdata, SEXP Nsnps, SEXP Nids, SEXP Coding, SEXP pedfilename,
		SEXP exportNumeric);

#ifdef __cplusplus
}
#endif


#endif

