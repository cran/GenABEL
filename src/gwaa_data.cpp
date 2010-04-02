#include <R.h>
#include <stdio.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP get_int_snp_matrix(SEXP Nids, SEXP Nsnps, SEXP data, SEXP Transposed)
	{
		int msk[4] = {192,48,12,3};
		int ofs[4] = {6,4,2,0};
		int i,j,m,idx=0;
		char str;
		int nsnps = INTEGER(Nsnps)[0];
		int nids = INTEGER(Nids)[0];
		int nbytes;
		bool transpose = (bool) LOGICAL(Transposed)[0];
		SEXP outdata;

		//		if ((nsnps % 4) == 0) nbytes = nsnps/4; else nbytes = ceil(1.*nsnps/4.);
		nbytes = unsigned(ceil(double(nids)/4.) + 0.5);

		if (transpose)
			PROTECT(outdata = allocMatrix(INTSXP, nsnps, nids));
		else
			PROTECT(outdata = allocMatrix(INTSXP, nids, nsnps));

		for (m=0;m<nsnps;m++) {
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = (char) RAW(data)[m*nbytes + i];
				//				Rprintf("a = %c\n",str);
				//				Rprintf("a = %u\n",(unsigned short int) str);
				for (j=0;j<4;j++) {
					int gt = str & msk[j];
					gt  >>= ofs[j];
					gt--;
					//Rprintf("%d (%d) ",idx,gt-1);
					if (transpose) {
						//Rprintf("m=%d nsnps=%d idx=%d (m+nsnps*idx)=%d gt=%d\n",m,nsnps,idx,m+nsnps*idx,gt);
						INTEGER(outdata)[m+nsnps*idx] = gt;
						if (gt<0) INTEGER(outdata)[m+idx*nsnps] = NA_INTEGER;
					} else {
						INTEGER(outdata)[m*nids+idx] = gt;
						if (gt<0) INTEGER(outdata)[m*nids+idx] = NA_INTEGER;
					}
					idx++;
					if (idx>=nids) {idx=0;break;}
				}
			}
		}

		UNPROTECT(1);

		return(outdata);
	}


	SEXP get_impute_snp_matrix(SEXP Nids, SEXP Nsnps, SEXP data) //, SEXP AlleleFreq)
	{
		int msk[4] = {192,48,12,3};
		int ofs[4] = {6,4,2,0};
		int i,j,m,idx=0;
		char str;
		int nsnps = INTEGER(Nsnps)[0];
		int nids = INTEGER(Nids)[0];
		int nbytes;
		SEXP outdata;

		//		if ((nsnps % 4) == 0) nbytes = nsnps/4; else nbytes = ceil(1.*nsnps/4.);
		nbytes = unsigned(ceil(double(nids)/4.) + 0.5);

		PROTECT(outdata = allocMatrix(REALSXP, nsnps, (nids*3)));

		for (m=0;m<nsnps;m++) {
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = (char) RAW(data)[m*nbytes + i];
				//				Rprintf("a = %c\n",str);
				//				Rprintf("a = %u\n",(unsigned short int) str);
				for (j=0;j<4;j++) {
					int gt = str & msk[j];
					gt  >>= ofs[j];
					gt--;
					//Rprintf("%d (%d) ",idx,gt-1);
					//Rprintf("m=%d nsnps=%d idx=%d (m+nsnps*idx)=%d gt=%d\n",m,nsnps,idx,m+nsnps*idx,gt);
					REAL(outdata)[m+nsnps*(idx*3)] = 0.;
					REAL(outdata)[m+nsnps*(idx*3)+nsnps*1] = 0.;
					REAL(outdata)[m+nsnps*(idx*3)+nsnps*2] = 0.;
					/**
					if (gt<0) {
						double q = REAL(AlleleFreq)[m];
						if (q<0.001) q=0.001;
						if (q>0.999) q=0.999;
						REAL(outdata)[m+nsnps*(idx*3)] = (1.-q)*(1.-q);
						REAL(outdata)[m+nsnps*(idx*3)+nsnps*1] = 2.*q*(1.-q);
						REAL(outdata)[m+nsnps*(idx*3)+nsnps*2] = q*q;
					} else {
					**/
						if (gt == 0) REAL(outdata)[m+nsnps*(idx*3)] = 1.;
						if (gt == 1) REAL(outdata)[m+nsnps*(idx*3)+nsnps*1] = 1.;
						if (gt == 2) REAL(outdata)[m+nsnps*(idx*3)+nsnps*2] = 1.;
//					}
					idx++;
					if (idx>=nids) {idx=0;break;}
				}
			}
		}

		UNPROTECT(1);

		return(outdata);
	}



#ifdef __cplusplus
}
#endif
