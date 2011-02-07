#include <R.h>
#include <Rinternals.h>



#ifdef __cplusplus
extern "C" {
#endif



void snp_summary_exhweWrapper(double *indata, unsigned long int indataHeight,
    unsigned long int indataWidth,	double *outdata,
	unsigned long int &outdataNcol, unsigned long int &outdataNrow,	unsigned int narg,
	SEXP *argList);

void snp_summary_exhwe_Processor(unsigned int *gt, unsigned int nids, double *out);


#ifdef __cplusplus
}
#endif

