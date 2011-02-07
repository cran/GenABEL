#ifndef __DAUTIL_H__
#define __DAUTIL_H__

#include <Rdefines.h>

#define error_R Rprintf("ERROR in Rstaff:"); Rprintf

#ifdef __cplusplus
extern "C" {
#endif

SEXP extract_text_file_column_cpp(SEXP, SEXP);

#ifdef __cplusplus
}
#endif


#endif
