#include "Rstaff.h"
#include <Rdefines.h>

#include "FilteredMatrix.h"
#include "Logger.h"

#ifdef __cplusplus
extern "C" {
#endif

    static void FilteredMatrixRFinalizer(SEXP x) {
   		if (x == R_NilValue) return;
   		FilteredMatrix* p = (FilteredMatrix *) EXTPTR_PTR(x);
   		if (p == NULL) return;
   		wrapperLog << "Finalizing FilteredMatrix: "<< (void *)p << endl;
   		delete p;
   	}

    static void FilteredAndAbstractRFinalizer(SEXP x) {
   		if (x == R_NilValue) return;
   		FilteredMatrix* p = (FilteredMatrix *) EXTPTR_PTR(x);
   		if (p == NULL) return;
   		AbstractMatrix *nestedMatrix = p->getNestedMatrix();
   		wrapperLog << "Finalizing FilteredMatrix: "<< (void *)p <<
   		" and nested AbstractMatrix "<< (void*)nestedMatrix << endl;
   		delete p;
   		delete nestedMatrix;            
   	}

   	SEXP setFilteredArea_R(SEXP filteredMatrixSEXP, SEXP selectedRows, SEXP selectedCols) {
        unsigned long i;

        vector<unsigned long> rowMask;
        for(i=0;i<((unsigned long) length(selectedRows));i++)
            rowMask.push_back(INTEGER(selectedRows)[i]-1);

        vector<unsigned long> colMask;
        for(i=0;i<((unsigned long) length(selectedCols));i++)
            colMask.push_back(INTEGER(selectedCols)[i]-1);

        checkPointer(filteredMatrixSEXP);
	    FilteredMatrix *filteredMatrix = (FilteredMatrix*)R_ExternalPtrAddr(filteredMatrixSEXP);
	    filteredMatrix->setFilteredArea(rowMask,colMask);

		return filteredMatrixSEXP;
    }

	SEXP create_FilteredMatrixFromAbstractMatrix_R(SEXP abstractMatrixSEXP) {

	    AbstractMatrix *abstractMatrix = (AbstractMatrix*)R_ExternalPtrAddr(abstractMatrixSEXP);

		FilteredMatrix * p = NULL;

		try {
			p = new FilteredMatrix(*abstractMatrix);
			cout << "create_FilteredMatrixFromAbstractMatrix_R = " << (void*)p << endl; 
		} catch (int errcode) {
			return R_NilValue;
		}

		if (p == NULL) {
			error_R("Error creating FilteredMatrix.\n");
			return R_NilValue;
		}

		SEXP val = R_MakeExternalPtr(p, install("FilteredMatrix"), R_NilValue);
		R_RegisterCFinalizerEx(val, FilteredMatrixRFinalizer, (Rboolean) TRUE);
		return val;
	}


	SEXP create_FilteredMatrixFromFilteredMatrix_R(SEXP filteredMatrixSEXP) {

	    FilteredMatrix *filteredMatrix = (FilteredMatrix*)R_ExternalPtrAddr(filteredMatrixSEXP);
		FilteredMatrix * p = NULL;

		try {
			p = new FilteredMatrix(*filteredMatrix);
			cout << "create_FilteredMatrixFromFilteredMatrix_R = " << (void*)p << endl; 
		} catch (int errcode) {
			return R_NilValue;
		}

		if (p == NULL) {
			error_R("Error creating FilteredMatrix.\n");
			return R_NilValue;
		}
		SEXP val = R_MakeExternalPtr(p, install("FilteredMatrix"), R_NilValue);
		R_RegisterCFinalizerEx(val, FilteredMatrixRFinalizer, (Rboolean) TRUE);
		return val;
	}

    SEXP disconnectFiltered_R(SEXP s) {
   		FilteredMatrixRFinalizer(s);
   		R_ClearExternalPtr(s);
   		return R_NilValue;
   	}

    SEXP disconnectFilteredAndAbstract_R(SEXP s) {
    	FilteredAndAbstractRFinalizer(s);
   		R_ClearExternalPtr(s);
   		return R_NilValue;
   	}

    SEXP open_FilteredMatrix_R(SEXP fname, SEXP cacheMb, SEXP ReadOnly) {
   		unsigned long cachesizeMb = (unsigned long) INTEGER(cacheMb)[0];
   		bool readonly = LOGICAL(ReadOnly)[0];
   		string filename = CHAR(STRING_ELT(fname,0));

   		if (cachesizeMb<0) {
   			error_R("negative cache size\n");
   			return R_NilValue;
   		}

   		FilteredMatrix* fm = NULL;

   		try {
   			FileVector *fv = new FileVector(filename,cachesizeMb,readonly);
   			fm = new FilteredMatrix(*fv);
//   			cout << "open_FilteredMatrix_R, ptr = " << (long)fm << endl;
   		} catch (int errcode) {
   			return R_NilValue;
   		}

   		if (fm == NULL) {
   			error_R("pointer is NULL\n");
   			return R_NilValue;
   		}
   		SEXP val = R_MakeExternalPtr(fm, Rf_install("FilteredMatrix"), R_NilValue);
   		R_RegisterCFinalizerEx(val, FilteredAndAbstractRFinalizer, (Rboolean) TRUE);
   		return val;
   	}


#ifdef __cplusplus
}
#endif
