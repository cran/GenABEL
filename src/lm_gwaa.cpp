/**
 *
 * 2008.05.05 by Yurii Aulchenko, EMCR
 *
 **/

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
using namespace std; 
#include <R.h>
#include "mematrix.h"
#include "mematri1.h"
#include "reg1data.h"
#include "reg1.h"
#include "mematriR.h"


void getgtvec(char *gtdata, int *gtint, int nids, int nbytes, int whichsnp)
{
  int msk[4] = {192,48,12,3};
  int ofs[4] = {6,4,2,0};
  int i,j,idx=0;
  char str;
  for (i=0;i<nbytes;i++) {
    str = gtdata[whichsnp*nbytes + i];
    for (j=0;j<4;j++) {
      gtint[idx] = str & msk[j]; 
      gtint[idx] >>= ofs[j];
      gtint[idx++]--;
      if (idx>=nids) {break;}
    }
  }
}

// gtmode: 1=additive, 2=dominant 3=recessive 4=overdominant
void convert_gt(int *gtint, int nids, int gtmode)
{
	if (gtmode==2) {
		for (int i=0;i<nids;i++) if (gtint[i]==2) gtint[i]=1;
	} else if (gtmode==3) {
		for (int i=0;i<nids;i++) if (gtint[i]==1) gtint[i]=0; else if (gtint[i]==2) gtint[i]=1;
	} else if (gtmode==4) {
		for (int i=0;i<nids;i++) if (gtint[i]==2) gtint[i]=0;
	}
}

extern "C" {
void linreg_gwaa(double *Y, double *X, char *gtdata, int *Nids, int *NXcol, int *Nsnps, int *GTMode, double *out)
   {
	int nids = (*Nids);
	int nsnps = (*Nsnps);
	int nxcol = (*NXcol);
	int gtmode = (*GTMode);
	int * igt = new (nothrow) int[nids];
	int nbytes;
//	if ((nsnps % 4) == 0) nbytes = nids/4; else nbytes = int(ceil(1.*double(nids)/4.));
	nbytes = int(ceil(1.*double(nids)/4.));
	for (int snp=0;snp<nsnps;snp++)
	{
		getgtvec(gtdata,igt,nids,nbytes,snp);
		if (gtmode>1) convert_gt(igt,nids,gtmode);
		regdata cdata(Y,X,igt,nids,nxcol,1);
		if (cdata.nids <= 1 || cdata.ismono) {
			out[snp] = cdata.nids;
			out[nsnps+snp] = -999.9;
			out[2*nsnps+snp] = -999.9;
		} else {
			linear_reg creg(cdata,0);
			out[snp] = cdata.nids;
			out[nsnps+snp] = creg.beta[(creg.beta).nrow-1];
			out[2*nsnps+snp] = creg.sebeta[(creg.beta).nrow-1];
		}
	}
	delete [] igt;
	return;
   }
}

extern "C" {
void logreg_gwaa(double *Y, double *X, char *gtdata, int *Nids, int *NXcol, int *Nsnps, int *GTMode, double *out)
   {
	int nids = (*Nids);
	int nsnps = (*Nsnps);
	int nxcol = (*NXcol);
	int gtmode = (*GTMode);
	int * igt = new (nothrow) int[nids];
	int nbytes;
//	if ((nsnps % 4) == 0) nbytes = nids/4; else nbytes = int(ceil(1.*double(nids)/4.));
	nbytes = int(ceil(1.*double(nids)/4.));
	for (int snp=0;snp<nsnps;snp++)
	{
		getgtvec(gtdata,igt,nids,nbytes,snp);
		if (gtmode>1) convert_gt(igt,nids,gtmode);
		regdata cdata(Y,X,igt,nids,nxcol,1);
		if (cdata.nids <= 1 || cdata.ismono) {
			out[snp] = cdata.nids;
			out[nsnps+snp] = -999.9;
			out[2*nsnps+snp] = -999.9;
		} else {
			logistic_reg creg(cdata,0,7);
			out[snp] = cdata.nids;
			out[nsnps+snp] = creg.beta[(creg.beta).nrow-1];
			out[2*nsnps+snp] = creg.sebeta[(creg.beta).nrow-1];
		}
	}
	delete [] igt;
	return;
   }
}

extern "C" {
void coxph_gwaa(double *Y, double *X, char *gtdata, int *Nids, int *NXcol, int *Nsnps, int *GTMode, double *out)
   {
	int nids = (*Nids);
	int nsnps = (*Nsnps);
	int nxcol = (*NXcol);
	int gtmode = (*GTMode);
	int * igt = new (nothrow) int[nids];
	int nbytes;
//	if ((nsnps % 4) == 0) nbytes = nids/4; else nbytes = int(ceil(1.*double(nids)/4.));
	nbytes = int(ceil(1.*double(nids)/4.));
	for (int snp=0;snp<nsnps;snp++)
	{
		getgtvec(gtdata,igt,nids,nbytes,snp);
		if (gtmode>1) convert_gt(igt,nids,gtmode);
		regdata cregdata(Y,X,igt,nids,nxcol,2);
		(cregdata.X).delete_column(0);
		if (cregdata.nids <= 1 || cregdata.ismono) {
			out[snp] = cregdata.nids;
			out[nsnps+snp] = -999.9;
			out[2*nsnps+snp] = -999.9;
		} else {
			coxph_data cdata(cregdata);
			coxph_reg creg(cdata,20,1.e-9,1.5e-12);
			out[snp] = cdata.nids;
			out[nsnps+snp] = creg.beta[(creg.beta).nrow-1];
			out[2*nsnps+snp] = creg.sebeta[(creg.beta).nrow-1];
		}
	}
	delete [] igt;
	return;
   }
}
