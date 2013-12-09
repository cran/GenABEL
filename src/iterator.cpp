#include <cstdarg>
#include "Rstuff.h"
#include "iterator_functions.h"
#include "iterator.h"
#include "gwaa_cpp.h"

#ifdef __cplusplus
extern "C" {
#endif

	MethodConvStruct methodConverter[] =
	{
			{ "sum", &sumWrapper },
			{ "prod", &prodWrapper },
			{ "sumpower", &sumpowerWrapper },
			{ "qtscore_glob", &qtscore_globWrapper },
			{ "snp_summary_exhwe", &snp_summary_exhweWrapper}
	};

	bool getDataReal(double *inData, unsigned int inDataRowLength, double *outData, unsigned int datasize,
			int step, unsigned long int index, unsigned int margin) {
		unsigned int k = 0;
		if (margin == 2) { // column-wise
			for (int i=0; i<step; i++) {
				for (unsigned long int j = 0; j < datasize; j++, k++) {
					//					cout << k << " " << index*datasize+k << " " << inData[index*datasize+k] << endl;
					outData[k] = inData[(index * inDataRowLength) + k];
				}
			}
		} else { // row-wise
			//Rprintf("Row %d ", index);
			for (int i=0; i<step; i++) {
				for (unsigned long int j = 0; j < datasize; j++, k++) {
					//					cout << k << " " << k*datasize+index << " " << inData[k*datasize+index] << endl;
					outData[k] = inData[(index + i) + (j * inDataRowLength)];
					//Rprintf("at %d: %f ", j, outData[j]);
				}
			}
			//Rprintf("\n\n");
		}
		return true;
	}

	bool getDataNew(AbstractMatrix *inData, double *outData, unsigned int datasize,
			int step,
			unsigned long int index, unsigned int margin) {
		if (margin == 2) { // column-wise
			try {
				for (int i=0;i<step;i++)
					inData->readVariableAs(index+i, &outData[i*datasize]);
			} catch (int errcode) {
				return false;
			}
		} else { // row-wise
			double dTmp;
			unsigned long int k=0;
			for (int i=0;i<step;i++) {
				for (unsigned long int j = 0; j < datasize; j++) {
					inData->readElementAs(j, index+i, dTmp);
					outData[k++] = dTmp;
				}
			}
		}
		return true;
	}

	void getDataOld(char const *inData, unsigned int inDataRowLength, double *outData,
			unsigned int datasize, int step,
			unsigned long int index, unsigned int margin) {
		int iTmp;
		char str;
		unsigned int k = 0;
		int msk[4] = { 192, 48, 12, 3 };
		int ofs[4] = { 6, 4, 2, 0 };
		unsigned short int last_till = 4;
		unsigned int nbytes; // always the length of a compressed row!
		if ((inDataRowLength % 4) == 0) {
			nbytes = (unsigned int) (inDataRowLength / 4);
		} else {
			nbytes = (unsigned int) ceil(1. * inDataRowLength / 4.);
			last_till = (unsigned short int) (inDataRowLength % 4);
		}

		double zero = 0.0;

		if (margin == 2) { // row-wise
			for (int iii = 0; iii < step; iii++) {
				unsigned int offset = (index + iii) * nbytes;
				for (unsigned int i = offset; i < offset + nbytes; i++) {
					str = inData[i];
					unsigned short int till = 4;
					if (i == offset + nbytes - 1) till = last_till;
					for (int j = 0; j < till; j++, k++) {
						iTmp = str & msk[j];
						iTmp >>= ofs[j];
						outData[k] = static_cast<double> (iTmp);
						outData[k] -= 1;
						if (iTmp == 0) outData[k] = 0 / zero;
					}
				}
			}
		} else { // column-wise
			int j = index % 4; // index inside char/byte
			//Rprintf("Index: %d, j: %d\n\n", index, j);
			unsigned int offset = (int)(index / 4); // starting point in char array
			unsigned int i;
			for (int iii = 0; iii < step; iii++) {
				for (unsigned int counter = 0; counter < datasize; counter++, k++) {
					i = (counter * nbytes) + offset;

					int byteindex, bitpairindex;
					if (j + iii <= 3) {
						byteindex = i;
						bitpairindex = j + iii;
					} else {
						int byteshift = (j + iii) / 4;
						int rest = (j + iii) % 4;
						byteindex = i + byteshift;
						bitpairindex = rest;
					}

					//Rprintf("At %d:%d indices: %d:%d ", i, j, byteindex, bitpairindex);
					str = inData[byteindex];
					iTmp = str & msk[bitpairindex];
					iTmp >>= ofs[bitpairindex];
					outData[k] = static_cast<double> (iTmp);
					outData[k] -= 1;
					if (iTmp == 0) outData[k] = 0 / zero;
					//Rprintf("At %d: %d ", k, (int)(outData[k]));
					//if (insloc > 0 && insloc % 10 == 0) Rprintf("\n");
				}
			}
			//Rprintf("\n");
			//Rprintf("Sumtotal at %d: %d\n", index, (int)sumtotal);
		}
	}


	SEXP iteratorGA(SEXP primedata, SEXP nrids, SEXP nrobs, SEXP method, SEXP outputtype,
			SEXP margin, SEXP nrstep, SEXP nrarg, ...) {

		// Check and get the data supplied
		unsigned long int nids, nobs;
		unsigned int intype = 0;
		// 0 = "databel", 1 = "snp.data", 2 = matrix
		AbstractMatrix *pDataNew = NULL;
		char const *pDataOld = NULL;
		double *pDataReal = NULL;
		int step = INTEGER(nrstep)[0];

		if (step<1) {
			error_R("nrstep < 1\n");
			return R_NilValue;
		}

		//Rprintf("TYPEOF(data)=%d\n",TYPEOF(primedata));
		SEXP Rclass = GET_CLASS(primedata);
		string klass;
		if (Rclass != R_NilValue)
			klass = CHAR(STRING_ELT(Rclass,0));
		else
			klass = "matrix";
		//Rprintf("GET_CLASS(data)=%s\n",klass.c_str());

		SEXP data;

		if (klass == "databel") {
			SEXP dataslot;
			PROTECT(dataslot = allocVector(STRSXP, (R_len_t) 1));
			SET_STRING_ELT(dataslot, 0, mkChar("data"));
			data = R_do_slot(primedata, dataslot);
			UNPROTECT(1);
			if (TYPEOF(data) != EXTPTRSXP) {
				error_R("Class of data is 'databel', but slot 'data' does not contain EXTPTRSXP\n");
				return R_NilValue;
			}
			intype = 0;
			pDataNew = getAbstractMatrixFromSEXP(data);
			if (pDataNew == NULL) {
				error_R("Pointer to data is NULL\n");
				return R_NilValue;
			}
			nids = pDataNew->getNumVariables();
			nobs = pDataNew->getNumObservations();
		} else if (klass == "snp.mx") {
			data = primedata;
			if (TYPEOF(data) != RAWSXP) {
				error_R("Class of data is 'snp.mx', but does not contain RAWSXP\n");
				return R_NilValue;
			}
			intype = 1;
			pDataOld = (char const *) RAW(data);
			//			nids = INTEGER(nrids)[0];
			//			nobs = INTEGER(nrobs)[0];
			nobs = INTEGER(nrids)[0];
			nids = INTEGER(nrobs)[0];
		} else if (klass == "matrix") {
			// DEALING WITH REGULAR MATRICES HERE
			data = primedata;
			if (TYPEOF(data) != REALSXP) {
				error_R("Class of data is 'matrix', but does not contain REALSXP\n");
				return R_NilValue;
			}
			intype = 2;
			pDataReal = (double *) REAL(data);
			nobs = INTEGER(nrids)[0];
			nids = INTEGER(nrobs)[0];
		} else {
			error_R("Incorrect data type\n");
			return R_NilValue;
		}

		// Find out and check the function supplied
		char const *methodName = CHAR(STRING_ELT(method, 0));
		myfunctiontype *pMethod = NULL;
		for (unsigned int i = 0; i < sizeof(methodConverter); i++) {
			if (strcmp(methodConverter[i].methodName, methodName) == 0) {
				pMethod = methodConverter[i].functionPtr;
				break;
			}
		}
		if (pMethod == NULL) {
			error_R("No (valid) function supplied\n");
			return R_NilValue;
		}

		// Find out the desired output type (file) supplied
		bool fv = true;
		char const *outputName = CHAR(STRING_ELT(outputtype, 0));
		if (strcmp(outputName, "R") == 0) {
			fv = false;
		}

		// Get the margin supplied
		int mar = INTEGER(margin)[0];
		if (mar < 1 || mar > 2) {
			error_R("No (valid) margin supplied\n");
			return R_NilValue;
		}

		// Get the nr. of additional arguments supplied
		unsigned int narg = INTEGER(nrarg)[0];

		// Get the additional parameters supplied, if any
		unsigned int argListSize = narg > 0 ? narg : 1;
		SEXP * argList = new (std::nothrow) SEXP [argListSize];
		if (argList == NULL) {error_R("can not allocate RAM for argList\n");return R_NilValue;}
		va_list ap;
		va_start(ap, nrarg); // nrarg is last known parameter
		for (unsigned int i = 0; i < narg; i++) {
			argList[i] = va_arg(ap, SEXP);
		}
		va_end(ap);

		// The actual data handling part:

		unsigned long int ncol, nrow;
		if (mar == 1) { // row-wise
			ncol = nobs;
			nrow = nids;
		} else { // column-wise (default)
			ncol = nids;
			nrow = nobs;
		}

		if ((ncol % step) != 0) {
			error_R("ncol not divisable by step\n");
			delete [] argList;
			return R_NilValue;
		}

		//cout << "before cycle 0:" << ncol << " " << nids << endl;
		unsigned long int nrow_new, ncol_multi;

		// Get the dimensions of the output the function of our choosing will be giving
		pMethod(0, nrow, step, 0, ncol_multi, nrow_new, narg, argList);
		//Rprintf("ncol_multi,nrow_new=%d,%d\n",ncol_multi,nrow_new);
		//cout << "ncol_multi=" << ncol_multi << endl;
		//cout << "nrow_new=" << nrow_new << endl;
		// Allocate vector
		// Start output SEXP for passing to R
		// Even when the output is put into a filevector, we still return an (empty) SEXP
		SEXP out;
		// Declare output filevector (whether we'll be using it or not)
		FileVector* tmpFV = NULL;
		FilteredMatrix* outFV = NULL;
		//AbstractMatrix * outFV;
		unsigned long int ncol_cs = (unsigned long int) ncol/step;
		if (!fv) {
			// Initialize output matrix once real number of rows is known
			// ASSUMPTION: nrow_new remains constant over calls to function wrapper
			// PROTECT(out = allocVector(REALSXP, (R_len_t)(ncol_cs * ncol_multi
			//		* nrow_new)));
			PROTECT(out = allocMatrix(REALSXP, (R_len_t)ncol_cs * nrow_new,
					(R_len_t)(ncol_multi)));
		} else {
			// To avoid null pointer error, make output SEXP to return (although it will be empty)
			PROTECT(out = allocVector(REALSXP, (R_len_t) 1));
			try {
				initializeEmptyFile(outputName, ncol_cs * nrow_new,
						ncol_multi, FLOAT, false);
			} catch (int errcode) {
				error_R("Failed in iterator - call - initializeEmptyFile");
				delete [] argList;
				UNPROTECT(1);
				return R_NilValue;
			}
			try {
				tmpFV = new FileVector(outputName,64,false);
				outFV = new FilteredMatrix(*tmpFV);
				//outFV = new FileVector(outputName, 64, false);
			} catch (int errcode) {
				error_R("Cannot initialize output file\n");
				delete [] argList;
				UNPROTECT(1);
				return R_NilValue;
			}
		}

		double * internal_data = new (std::nothrow) double [nrow*step];
		if (internal_data == NULL) {
			error_R("can not allocate RAM for internal_data\n");
			delete [] argList;
			UNPROTECT(1);
			return R_NilValue;
		}
		double * out_data = new (std::nothrow) double [nrow_new * ncol_multi];
		//Rprintf("ncol_multi,nrow_new=%d,%d\n",ncol_multi,nrow_new);
		if (out_data == NULL) {
			error_R("can not allocate RAM for out_data\n");
			delete [] argList;
			delete [] internal_data;
			UNPROTECT(1);
			return R_NilValue;
		}

		// Read in data and apply function (row- or column-wise)
		// use 'STEP' here?
		//cout << "before cycle:" << ncol << " " << step << endl;

		int rowlength = nrow;
		if (mar == 1) {
			rowlength = ncol;
		}

		for (unsigned long int i = 0; i < ncol; i+=step) {
			//cout << "in cycle: " << i << endl;
			// Get row or column
			if (intype==0) {
				getDataNew(pDataNew, internal_data, nrow, step, i, mar);
			} else if (intype==1) {
				getDataOld(pDataOld, rowlength, internal_data, nrow, step, i, mar);
			} else if (intype==2) {
				getDataReal(pDataReal, rowlength, internal_data, nrow, step, i, mar);
			} else {
				error_R("unsupported intype\n");
				delete [] argList;
				delete [] internal_data;
				UNPROTECT(1);
				return R_NilValue;
			}

			/**
			cout << "internal_data " << i << ":";
			for (int iii=0;iii<step;iii++)
				for (int jjj=0;jjj<nrow;jjj++) cout << " " << internal_data[iii*nrow+jjj];
			cout << endl;
			 **/

			// Apply function of choosing
			pMethod(internal_data, nrow, step, out_data, ncol_multi, nrow_new, narg,
					argList);
			//cout << "pMethod returns " << out_data[0] << endl;

			// Write analyzed data to R vector or filevector
			unsigned long int i_cs = (unsigned long int) (i / step);
			for (unsigned long int j = 0; j < ncol_multi; j++) {
				if (!fv) {
					for (unsigned long int k = 0; k < nrow_new; k++) {
						// Add to output SEXP
						//cout << "set element " << (i_cs * ncol_multi + j) * nrow_new + k <<
						//		"from k " << k << " to value " << out_data[k] << endl;
						//REAL(out)[(i_cs * ncol_multi + j) * nrow_new + k]
						//          = out_data[k];
						//REAL(out)[(i_cs * ncol_multi + j) * nrow_new + k]
						 //         = out_data[j*nrow_new+k];
						REAL(out)[j*ncol_cs+i_cs+k]
						          = out_data[j*nrow_new+k];
		}
				} else {
					outFV->writeVariableAs(i_cs * ncol_multi + j, &out_data[j*nrow_new]);
				}
			}
		}

		if (fv) {
			//tmpFV.deInitialize();
			delete outFV;
			delete tmpFV;
		}

		delete [] argList;
		delete [] internal_data;
		delete [] out_data;

		UNPROTECT(1);
		return out;
	}


	/**
	// OLD STUFF BELOW HERE:
	// iterator and other stuff
	SEXP databel_impute_prob_2_databel_mach_dose(SEXP imputedata, SEXP OutFileName, SEXP CacheSizeMb)
	{
		CHECK_PTR(imputedata);
		//Rprintf("CHECK PASSED\n");
		AbstractMatrix * p = (AbstractMatrix*)R_ExternalPtrAddr(imputedata);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}
		std::string outfilename = CHAR(STRING_ELT(OutFileName,0));
		unsigned long int cachesizeMb = (unsigned long int) INTEGER(CacheSizeMb)[0];
		if (cachesizeMb<0) {
			error_R("negative cache size\n");
			return R_NilValue;
		}

		unsigned int nvars, nobss;

		try {
			nvars = (unsigned int) p->getNumVariables();
		} catch (int errcode) {
			error_R("ERROR: can not getNumVariables\n");
			return R_NilValue;
		}
		try {
			nobss = (unsigned int) p->getNumObservations();
		} catch (int errcode) {
			error_R("ERROR: can not getNumObservations\n");
			return R_NilValue;
		}

		if ((nobss % 3) != 0)
		{
			error_R("ERROR: nobss not divisable by 3\n");
			return R_NilValue;
		}

		unsigned int new_nvars = nvars, new_nobss = (unsigned int) nobss/3;

		try {
			// last flag -- override
			initializeEmptyFile(outfilename, new_nvars, new_nobss, FLOAT, false);
		} catch (int errcode) {
			error_R("failed in databel_impute_prob_2_databel_mach_dose - call - initializeEmptyFile");
			return R_NilValue;
		}

		Rprintf("READY TO INI NEW\n");
		AbstractMatrix * outFV;
		try {
			outFV = new FileVector(outfilename,cachesizeMb);
		} catch (int errcode) {
			error_R("can not ini outfile\n");
			return R_NilValue;
		}

		float tmpout[new_nobss], tmpin[nobss];
		for (unsigned int var=0; var<nvars; var++)
		{
			p->readVariableAs(var, tmpin);
			//for (int k=0;k<nobss;k++) Rprintf("%f ",tmpin[k]); Rprintf("\n");
			unsigned int j = 0;
			for (unsigned int obs=0;obs<nobss;obs+=3)
			{
				tmpout[j++] = 2.*tmpin[obs+2]+tmpin[obs+1];
			}
			//for (int k=0;k<new_nobss;k++) Rprintf("%f ",tmpout[k]); Rprintf("\n");
			outFV->writeVariableAs(var,tmpout);
		}

		delete outFV;
		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		UNPROTECT(1);
		return ret;
	}


	SEXP databel_impute_prob_2_databel_mach_prob(SEXP imputedata, SEXP OutFileName, SEXP CacheSizeMb)
	{
		CHECK_PTR(imputedata);
		AbstractMatrix * p = (AbstractMatrix*)R_ExternalPtrAddr(imputedata);
		if (p == NULL) {
			error_R("pointer is NULL\n");
			return R_NilValue;
		}
		std::string outfilename = CHAR(STRING_ELT(OutFileName,0));
		unsigned long int cachesizeMb = (unsigned long int) INTEGER(CacheSizeMb)[0];
		if (cachesizeMb<0) {
			error_R("negative cache size\n");
			return R_NilValue;
		}

		unsigned int nvars, nobss;

		try {
			nvars = (unsigned int) p->getNumVariables();
		} catch (int errcode) {
			error_R("ERROR: can not getNumVariables\n");
			return R_NilValue;
		}
		try {
			nobss = (unsigned int) p->getNumObservations();
		} catch (int errcode) {
			error_R("ERROR: can not getNumObservations\n");
			return R_NilValue;
		}

		if ((nobss % 3) != 0)
		{
			error_R("ERROR: nobss not divisable by 3\n");
			return R_NilValue;
		}

		unsigned int new_nvars = 2*nvars, new_nobss = (unsigned int) nobss/3;

		try {
			// last flag -- override
			initializeEmptyFile(outfilename, new_nvars, new_nobss, FLOAT, false);
		} catch (int errcode) {
			error_R("failed in databel_impute_prob_2_databel_mach_prob - call - initializeEmptyFile");
			return R_NilValue;
		}

		AbstractMatrix * outFV;
		try {
			outFV = new FileVector(outfilename,cachesizeMb);
		} catch (int errcode) {
			error_R("can not ini outfile\n");
			return R_NilValue;
		}

		float tmpout1[new_nobss], tmpout2[new_nobss], tmpin[nobss];
		unsigned int coutvar = 0;
		for (unsigned int var=0; var<nvars; var++)
		{
			p->readVariableAs(var, tmpin);
			//for (int k=0;k<nobss;k++) Rprintf("%f ",tmpin[k]); Rprintf("\n");
			unsigned int j = 0;
			for (unsigned int obs=0;obs<nobss;obs+=3)
			{
				tmpout1[j] = tmpin[obs+1];
				tmpout2[j++] = tmpin[obs+2];
			}
			//for (int k=0;k<new_nobss;k++) Rprintf("%f ",tmpout2[k]); Rprintf("\n");
			//for (int k=0;k<new_nobss;k++) Rprintf("%f ",tmpout1[k]); Rprintf("\n");
			outFV->writeVariableAs(coutvar++,tmpout2);
			outFV->writeVariableAs(coutvar++,tmpout1);
		}

		delete outFV;
		SEXP ret;
		PROTECT(ret = allocVector(LGLSXP, 1));
		LOGICAL(ret)[0] = TRUE;
		UNPROTECT(1);
		return ret;
	}


	 **/

#ifdef __cplusplus
}
#endif
