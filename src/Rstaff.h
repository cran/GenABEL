#ifndef __RSTAFF_H__
#define __RSTAFF_H__

#include <string>
#include <cstring>
#include <R.h>


#include "convert_util.h"
#include "AbstractMatrix.h"
#include "FilteredMatrix.h"
#include "Transposer.h"
#include "dautil.h"


// maximal number of file-matrices allowed
// #define MAX_FM_OBJECTS 10

#ifdef __cplusplus
extern "C" {
#endif

	AbstractMatrix *getAbstractMatrixFromSEXP(SEXP s);

	//check if ptr valid
    void checkPointer(SEXP s);

#ifdef __cplusplus
}
#endif

#endif
