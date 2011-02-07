//#include "Rstaff.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "dautil.h"

#ifdef __cplusplus
extern "C" {
#endif

	SEXP extract_text_file_column_cpp(SEXP Filename, SEXP Whichcol)
	{
		std::string filename = CHAR(STRING_ELT(Filename,0));
		int ncols = (int) LENGTH(Whichcol);
		std::vector<int> whichcol(ncols);
		int maxcol = 0;
		for (int i=0;i<ncols;i++)
		{
			whichcol[i] = (int) INTEGER(Whichcol)[i];
			if (whichcol[i]>maxcol) maxcol=whichcol[i];
		}
		std::vector<std::string> outstrings;

		std::ifstream filetoread (filename.c_str());
		if (!filetoread) {
			error_R("can not open file '%s'\n\n",filename.c_str());
			return R_NilValue;
		}

		std::string tmpstr;
		while (getline(filetoread,tmpstr))
		{
			std::istringstream datas (tmpstr);
			std::string substr;
			std::vector<std::string> tmpvec;
			for (int i=0;i<=maxcol;i++) {
				datas >> substr;
				tmpvec.push_back(substr);
			}
			for (int i=0;i<ncols;i++) {
				//Rprintf("%i %i %s\n",i,whichcol[i],tmpvec[whichcol[i]].c_str());
				outstrings.push_back(tmpvec[whichcol[i]]);
			}
		}

		SEXP ret;
		PROTECT(ret = allocVector(STRSXP, (R_len_t) outstrings.size()));
		for (unsigned long int i = 0;i<outstrings.size();i++) SET_STRING_ELT(ret, i, mkChar(outstrings[i].c_str()));
		UNPROTECT(1);
		return ret;
	}


#ifdef __cplusplus
}
#endif
