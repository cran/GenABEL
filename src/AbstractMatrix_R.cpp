#include "Rstaff.h"

#include "FilteredMatrix.h"
#include "Logger.h"

// must be included after c++ headers!
#include <stdio.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif
	void checkPointer(SEXP s) {
		if (TYPEOF(s) != EXTPTRSXP) {
			errorLog << "Pointer is not EXTPTRSXP" << endl << errorExit;
		}
		if (R_ExternalPtrTag(s) != install("AbstractMatrix") && R_ExternalPtrTag(s) != install("FilteredMatrix")) {
			errorLog << "R_ExternalPtrTag(s) = " << (void*)R_ExternalPtrTag(s) << endl;
			errorLog << "Pointer is not AbstractMatrix nor FilteredMatrix" << endl << errorExit;
		}
	}

	AbstractMatrix *getAbstractMatrixFromSEXP(SEXP s){
		checkPointer(s);
		if (TYPEOF(s) == EXTPTRSXP) {
			return  ((AbstractMatrix*)R_ExternalPtrAddr(s))->castToAbstractMatrix();
		}
		errorLog << "External pointer not valid!" << endl << errorExit ;
		return NULL;
	}

	SEXP get_nvars_R(SEXP s) {
		//	    cout << "get_nvars_R()" << endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		SEXP out;
		PROTECT(out = allocVector(INTSXP, 1));
		unsigned int nvars = 0;

		try {
			nvars = (unsigned int) p->getNumVariables();
		} catch (int errcode) {
			nvars = 0;
		}

		if (nvars<=0) {
			out = R_NilValue;
		} else {
			INTEGER(out)[0] = nvars;
		}
		UNPROTECT(1);
		return out;
	}

	SEXP get_nobs_R(SEXP s) {
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		SEXP out;
		PROTECT(out = allocVector(INTSXP, 1));
		unsigned int nobss = 0;

		try {
			nobss = (unsigned int) p->getNumObservations();
		} catch (int errcode) {
			nobss = 0;
		}

		if (nobss<=0) {
			out = R_NilValue;
		} else {
			INTEGER(out)[0] = nobss;
		}
		UNPROTECT(1);
		return out;
	}

	SEXP setReadOnly_R(SEXP s, SEXP readOnly) {
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		bool readonly = LOGICAL(readOnly)[0];

		bool result = p->setReadOnly(readonly);

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = result?TRUE:FALSE;
		UNPROTECT(1);
		return ret;
	}

	SEXP get_all_varnames_R(SEXP s) {
		//	    testDbg << "get_all_varnames_R" << endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		//R_len_t nvars = (R_len_t) 0;
		unsigned long int nvars = 0;

		try {
			nvars = p->getNumVariables();
		} catch (int errcode) {
			return R_NilValue;
		}

		FixedChar tmp;
		SEXP ret;
		//cout << "get_all_varnames.nvars=" << nvars << endl;
		PROTECT(ret = allocVector(STRSXP, (R_len_t) nvars));
		//cout << "alloc done" << endl;

		try {
			for (unsigned long i = 0; i< nvars; i++) {
				tmp = p->readVariableName(i);
				SET_STRING_ELT(ret, i, mkChar(tmp.name));
			}
		} catch (int errcode) {
			error_R("something went terribly wrong in get_all_varnames_R\n");
			UNPROTECT(1);
			return ret;
		}
		UNPROTECT(1);
		return ret;
	}

	// !!!
	SEXP set_all_varnames_R(SEXP s, SEXP names) {
		//   	    testDbg << "set_all_varnames_R"<<endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

//		R_len_t nvars = (R_len_t) 0;
		unsigned long nvars = 0;

		try {
			nvars = p->getNumVariables();
		} catch (int errcode) {
			error_R("can not p->getNumVariables()\n");
			return R_NilValue;
		}

		// check that length of SEXP names is the same!!!

		for (unsigned long i = 0; i < nvars; i++) {
			string varname = CHAR(STRING_ELT(names,i));
			try {
				p->writeVariableName(i,FixedChar(varname));
			} catch (int errcode) {
				error_R("can not set variable name for variable %ul\n",i);
				return R_NilValue;
			}
		}

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		UNPROTECT(1);
		return ret;

	}

	SEXP get_all_obsnames_R(SEXP s) {
		//testDbg << "get_all_obsnames_R"<<endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		//R_len_t nobss = (R_len_t) 0;
		unsigned long int nobss = 0;

		try {
			nobss = p->getNumObservations();
		} catch (int errcode) {
			return R_NilValue;
		}

		FixedChar tmp;
		SEXP ret;
		PROTECT(ret = allocVector(STRSXP, (R_len_t) nobss));

		try {
			for (unsigned long i = 0; i< nobss; i++) {
				tmp = p->readObservationName(i);
				SET_STRING_ELT(ret, i, mkChar(tmp.name));
			}
		} catch (int errcode) {
			error_R("something went terribly wrong in get_all_obsnames_R\n");
			UNPROTECT(1);
			return ret;
		}
		UNPROTECT(1);
		return ret;
	}


	SEXP set_all_obsnames_R(SEXP s, SEXP names) {
		//testDbg << "set_all_obsnames_R"<<endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		//R_len_t nobss = (R_len_t) 0;
		unsigned long int nobss = 0;

		try {
			nobss = p->getNumObservations();
		} catch (int errcode) {
			error_R("can not p->getNumObservations()\n");
			return R_NilValue;
		}

		// check that length of SEXP names is the same!!!

		for (unsigned long i = 0; i < nobss; i++) {
			string obsname = CHAR(STRING_ELT(names,i));
			try {
				p->writeObservationName(i,FixedChar(obsname));
			} catch (int errcode) {
				error_R("can not set observation name for observation %ul\n",i);
				return R_NilValue;
			}
		}

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		UNPROTECT(1);
		return ret;

	}

	static void AbstractMatrixRFinalizer(SEXP x) {
		if (x == R_NilValue) return;
		AbstractMatrix* p = (AbstractMatrix *) EXTPTR_PTR(x);
		if (p == NULL) return;
		wrapperLog << "finalizing AbstractMatrix: " << (void*)p << endl;  
		delete p;
	}

	SEXP disconnect_R(SEXP s) {
		AbstractMatrixRFinalizer(s);
		R_ClearExternalPtr(s);
		return R_NilValue;
	}

	SEXP externalptr_is_null(SEXP s) {
		checkPointer(s);
		AbstractMatrix * p = (AbstractMatrix*)R_ExternalPtrAddr(s);
		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = FALSE;
		if (p == NULL) LOGICAL(ret)[0] = TRUE;
		UNPROTECT(1);
		return ret;
	}

	SEXP open_FileMatrix_R(SEXP fname, SEXP cacheMb, SEXP ReadOnly) {
		unsigned long cachesizeMb = (unsigned long) INTEGER(cacheMb)[0];
		bool readonly = LOGICAL(ReadOnly)[0];
		string filename = CHAR(STRING_ELT(fname,0));
		if (cachesizeMb<0) {
			error_R("negative cache size\n");
			return R_NilValue;
		}

		AbstractMatrix* p = NULL;

		try {
			p = new FileVector(filename,cachesizeMb,readonly);
			//cout << "open_FileMatrix_R, ptr = " << (void*)p << endl;
		} catch (int errcode) {
			return R_NilValue;
		}

		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}
		SEXP val = R_MakeExternalPtr(p, Rf_install("AbstractMatrix"), R_NilValue);
		R_RegisterCFinalizerEx(val, AbstractMatrixRFinalizer, (Rboolean) TRUE);
		return val;
	}

	SEXP read_variable_double_FileMatrix_R(SEXP nvar, SEXP s) {
		//testDbg << "read_variable_float_FileMatrix_R"<<endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}
		unsigned long nvariable = (unsigned long) INTEGER(nvar)[0] - 1;
		unsigned int nobs = 0;
		try {
			nobs = p->getNumObservations();
		} catch (int errcode) {
			return R_NilValue;
		}
		double * internal_data = new (std::nothrow) double [nobs];

		try {
			p->readVariableAs(nvariable, internal_data);
		} catch (int errcode) {
			return R_NilValue;
		}

		SEXP out;
		PROTECT(out = allocVector(REALSXP, (R_len_t) p->getNumObservations()));
		for (unsigned long i=0;i< nobs; i++) REAL(out)[i] = internal_data[i];
		delete [] internal_data;

		UNPROTECT(1);

		return out;
	}

	/// magic number 9007199254740993, minimum value when (float)(x+1) != x

	SEXP write_variable_double_FileMatrix_R(SEXP nvar, SEXP data, SEXP s) {
		//testDbg << "write_variable_double_FileMatrix_R"<<endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}
		unsigned long nvariable = (unsigned long) INTEGER(nvar)[0] - 1;
		// here generally should be very careful -- what type of data is IN?

		unsigned int nvars = 0;
		unsigned int nobss = 0;

		try {
			nvars = p->getNumVariables();
		} catch (int errocode) {
			return R_NilValue;
		}

		if (nvariable <0 || nvariable >= nvars) {
			error_R("nvar (%lu) out of range!\n",nvariable);
			return R_NilValue;
		}

		try {
			nobss = p->getNumObservations();
		} catch (int errcode) {
			return R_NilValue;
		}

		double * internal_data = new (std::nothrow) double [nobss];

		if (internal_data == NULL) {
			error_R("internal_data pointer is NULL\n");
			return R_NilValue;
		}

		for (unsigned long i=0;i< nobss;i++) {
			internal_data[i] = REAL(data)[i];
		}

		try {
			p->writeVariableAs(nvariable, internal_data);
		} catch (int errcode) {
			delete [] internal_data;
			error_R("can not write variable %ul\n",nvariable);
			return R_NilValue;
		}

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		delete [] internal_data;

		UNPROTECT(1);
		return ret;
	}

	// !!!
	SEXP set_cachesizeMb_R(SEXP s, SEXP SizeMB)
	{
		//testDbg << "set_cachesizeMb_R"<<endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}
		unsigned long sizeMb = (unsigned long) INTEGER(SizeMB)[0];
		try {
			p->setCacheSizeInMb( sizeMb );
		} catch (int errcode) {
			error_R("cannot reset cache size\n");
			return R_NilValue;
		}

		SEXP ret;PROTECT(ret = allocVector(LGLSXP, 1));LOGICAL(ret)[0] = TRUE;UNPROTECT(1);
		return ret;

	}

	SEXP get_cachesizeMb_R(SEXP s)
	{
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		unsigned long sizeMb = 0;

		try {
			sizeMb = p->getCacheSizeInMb();
		} catch (int errcode) {
			return R_NilValue;
		}

		SEXP out;
		PROTECT(out = allocVector(INTSXP, 1));
		INTEGER(out)[0] = (int) sizeMb;
		UNPROTECT(1);
		return(out);
	}


	SEXP text2fvf_R(SEXP Fnames, SEXP IntPars)
	{
		string program_name = "text2fvf_R";
		string infilename = CHAR(STRING_ELT(Fnames,0));
		string outfilename = CHAR(STRING_ELT(Fnames,1));
		string rownamesfilename = CHAR(STRING_ELT(Fnames,2));
		string colnamesfilename = CHAR(STRING_ELT(Fnames,3));
		unsigned long rownames = (unsigned long) INTEGER(IntPars)[0];
		unsigned long colnames = (unsigned long) INTEGER(IntPars)[1];
		unsigned long skiprows = (unsigned long) INTEGER(IntPars)[2];
		unsigned long skipcols = (unsigned long) INTEGER(IntPars)[3];
		int transpose = (int) INTEGER(IntPars)[4];
		int Rmatrix = (int) INTEGER(IntPars)[5];
		unsigned short Type = (unsigned short ) INTEGER(IntPars)[6];
		string nanString = CHAR(STRING_ELT(Fnames,4));

		try {
			text2fvf(program_name,
					infilename, outfilename,
					rownamesfilename, colnamesfilename,
					rownames, colnames,
					skiprows, skipcols,
					transpose, Rmatrix, Type, false, nanString);
		} catch (int x) {
			error_R("failed in text2fvf_R\n");
			return R_NilValue;
		}

		SEXP ret;PROTECT(ret = allocVector(LGLSXP, 1));LOGICAL(ret)[0] = TRUE;UNPROTECT(1);
		return ret;
	}

	SEXP ini_empty_FileMatrix_R(SEXP fname, SEXP nvars, SEXP nobs, SEXP Type)
	{
		unsigned long numVariables = (unsigned long) INTEGER(nvars)[0];
		unsigned long nobservations = (unsigned long) INTEGER(nobs)[0];
		string filename = CHAR(STRING_ELT(fname,0));
		unsigned short int type = (unsigned short int) INTEGER(Type)[0];

		if (type <=0 || type > 8) {
			error_R("Unknown data type %u\n",type);
			return R_NilValue;
		}
		try {
			// last flag -- override
			initializeEmptyFile(filename, numVariables, nobservations, type, false);
		} catch (int errcode) {
			error_R("failed in ini_empty_FileMatrix_R");
			return R_NilValue;
		}

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		UNPROTECT(1);
		return ret;

	}

	//virtual void save(string newFilename, unsigned long nvars, unsigned long nobss, unsigned long * varindexes, unsigned long * obsindexes)
	SEXP save_R(SEXP New_file_name, SEXP IntPars, SEXP s)
	{
		//   dbg<<"save_R"<<endl;
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		string newFilename = CHAR(STRING_ELT(New_file_name,0));
		unsigned long nvars = (unsigned long) INTEGER(IntPars)[0];
		unsigned long nobss = (unsigned long) INTEGER(IntPars)[1];
		unsigned long * varindexes = new (std::nothrow) unsigned long [nvars];
		if (varindexes == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}
		unsigned long * obsindexes = new (std::nothrow) unsigned long [nobss];
		if (obsindexes == NULL) {
			error_R("pointer is NULL\n");
			delete [] varindexes;
			return R_NilValue;
		}

		for (unsigned long i = 0; i < nvars; i++)
			varindexes[i] = (unsigned long) INTEGER(IntPars)[i+2];
		for (unsigned long i = 0; i < nobss; i++) {
			obsindexes[i] = (unsigned long) INTEGER(IntPars)[i+2+nvars];
		}

		try {
			p->saveAs(newFilename,nvars,nobss,varindexes,obsindexes);
		} catch (int errcode) {
			error_R("can not save data to file %s\n",newFilename.c_str());
			delete [] obsindexes;
			delete [] varindexes;
			return R_NilValue;
		}

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		delete [] obsindexes;
		delete [] varindexes;

		UNPROTECT(1);
		return ret;
	}

	SEXP saveAsText(SEXP s, SEXP New_file_name, SEXP IntPars, SEXP NANString ) 	{
		AbstractMatrix * p = getAbstractMatrixFromSEXP(s);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}

		string newFilename = CHAR(STRING_ELT(New_file_name,0));
		string nanString = CHAR(STRING_ELT(NANString,0));
		bool showVarNames = LOGICAL(IntPars)[0];
		bool showObsNames = LOGICAL(IntPars)[1];
		bool transpose = LOGICAL(IntPars)[2];

        AbstractMatrix *transposed = p;
        string tmpFileName,tmpFileName2;
        if (!transpose){
		    Transposer transposer;
		    tmpFileName= p->getFileName() + string("_saveAsText_tmp");
		    tmpFileName2= p->getFileName() + string("_saveAsText_tmp2");
		    p->saveAs(tmpFileName);
		    transposer.process(tmpFileName, tmpFileName2, true);
		    transposed = new FileVector(tmpFileName2, p->getCacheSizeInMb());
		}

		try {
			transposed->saveAsText(newFilename, showVarNames, showObsNames, nanString);
		} catch (int errcode) {
			error_R("can not save data to file %s\n",newFilename.c_str());
			return R_NilValue;
		}
        if (!transpose){
            delete transposed;
            unlink(tmpFileName.c_str());
            unlink(tmpFileName2.c_str());
        }

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;

		UNPROTECT(1);
		return ret;
	}

	SEXP checkNumBits(){
	    if (sizeof(unsigned long) != 8) {
    		errorLog << "YOU APPEAR TO WORK ON 32-BIT SYSTEM. LARGE FILES ARE NOT SUPPORTED."<<endl;
    	}
	    return R_NilValue;
	}

	/**
	*  direction: 0 -- copy values from @values to @ptr,
	*             1 -- copy values from @ptr to @values...
	*/
	SEXP assignDoubleMatrix(SEXP ptr, SEXP obsIndexes, SEXP varIndexes, SEXP values, SEXP direction){
		//flush(cout);
		
	    unsigned long varNum, obsNum, obsIndexNum, varIndexNum;

		AbstractMatrix * p = getAbstractMatrixFromSEXP(ptr);
        double coeff = 1. * length(obsIndexes) / p->getNumObservations();

		unsigned long dir = (unsigned long) INTEGER(direction)[0];

        double *currentValues = 0;
        if (!(coeff < WRITE_SPEED_PROPORTION)) {
            currentValues = new double[p->getNumObservations()];
        }

        unsigned long varIndexesLength = length(varIndexes);
        unsigned long obsIndexesLength = length(obsIndexes);

	    for(varIndexNum = 0; varIndexNum < varIndexesLength; varIndexNum ++){
	        varNum = (unsigned long)INTEGER(varIndexes)[varIndexNum]-1;

	        if ( coeff < WRITE_SPEED_PROPORTION) {
		
    	        for(obsIndexNum = 0; obsIndexNum < obsIndexesLength; obsIndexNum ++){
	                obsNum = (unsigned long)INTEGER(obsIndexes)[obsIndexNum]-1;
	                try {
	                    if (dir==0) {
     	                    double value = REAL(values)[varIndexNum *obsIndexesLength + obsIndexNum];
                            p->writeElementAs(varNum, obsNum,value);
                        } else {
     	                    double value;
                            p->readElementAs(varNum, obsNum,value);
                            REAL(values)[varIndexNum *obsIndexesLength + obsIndexNum] = value;
                        }
                    } catch(int errorCode) {
                        return R_NilValue;
                    }
	            }
	        } else {
	            try {
                    if (dir==0) {
    	                p->readVariableAs(varNum, currentValues);
        	            for(obsIndexNum = 0; obsIndexNum < obsIndexesLength; obsIndexNum ++){
    	                    obsNum = (unsigned long)INTEGER(obsIndexes)[obsIndexNum] - 1;
        	                currentValues[obsNum] = REAL(values)[varIndexNum*obsIndexesLength+obsIndexNum];
        	            }
                        p->writeVariableAs(varNum, currentValues);
        	        } else {
    	                p->readVariableAs(varNum, currentValues);
        	            for(obsIndexNum = 0; obsIndexNum < obsIndexesLength; obsIndexNum ++){
    	                    obsNum = (unsigned long)INTEGER(obsIndexes)[obsIndexNum] - 1;
        	                REAL(values)[varIndexNum*obsIndexesLength+obsIndexNum] = currentValues[obsNum];
        	            }
        	        }
                } catch(int errorCode){
                    delete [] currentValues;
                    return R_NilValue;
                }
	        }
	    }

        if (!(coeff < WRITE_SPEED_PROPORTION)) {
	        delete [] currentValues;
	    }

		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		UNPROTECT(1);
		//flush(cout);
		return ret;
	}



#ifdef __cplusplus
}
#endif
//
// OLD STRANGE STAFF
//
#ifdef __cplusplus
extern "C" {
#endif

	//        .Fortran("dqrls",
	//                  qr = x, n = n, p = p,
	//                  y = tra, ny = ny,
	//                  tol = as.double(tol),
	//                  coefficients = mat.or.vec(p, ny),
	//                  residuals = y, effects = y, rank = integer(1L),
	//                  pivot = 1L:p, qraux = double(p), work = double(2*p),
	//                  PACKAGE="base")$coefficients[2]

	void dqrls_(
			double*, int*, int*,
			double*, int*,
			double*,
			double*,
			double*, double*, int*,
			int*, double*, double*
	);

	//    .Fortran("ch2inv", x = x, nr, size, v = matrix(0, nrow = size,
	//        ncol = size), info = integer(1L), DUP = FALSE, PACKAGE = "base")
	void ch2inv_(
			double*,int*,int*,double*,int*
	);


	//      subroutine dqrls(x,n,p,y,ny,tol,b,rsd,qty,k,jpvt,qraux,work)
	//      integer n,p,ny,k,jpvt(p)
	//      double precision x(n,p),y(n,ny),tol,b(p,ny),rsd(n,ny),
	//     .                 qty(n,ny),qraux(p),work(p)

	void CPP_dqrls(
			double * x, int * n, int * p,
			double * y, int * ny,
			double * tol,
			double * b,
			double * rsd, double * qty, int * k,
			int * jpvt, double * qraux, double * work
	)
	{
		dqrls_(x,n,p,y,ny,tol,b,rsd,qty,k,jpvt,qraux,work);
	}

#ifdef __cplusplus
}
#endif

