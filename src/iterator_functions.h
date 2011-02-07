#ifndef __ITERATOR_FUNCTIONS_H__
#define __ITERATOR_FUNCTIONS_H__

#ifdef __cplusplus
extern "C" {
#endif

	typedef void (myfunctiontype)(double *, unsigned long int, unsigned long int,
			double *, unsigned long int &, unsigned long int &, unsigned int, SEXP *);

	struct MethodConvStruct {
		const char *methodName;
		myfunctiontype *functionPtr;
	};

	double prod(double *, unsigned int );
	void prodWrapper(double *, unsigned long int , unsigned long int , double *,
					unsigned long int &, unsigned long int &,
					unsigned int , SEXP *);
	double sum(double *, unsigned int , bool );
	void sumWrapper(double *, unsigned long int , unsigned long int , double *,
					unsigned long int &, unsigned long int &,
					unsigned int , SEXP *);
	double sumpower(double *, unsigned int , int );
	void sumpowerWrapper(double *, unsigned long int , unsigned long int ,
					double *, unsigned long int &,
					unsigned long int &, unsigned int , SEXP *);

	void qtscore_globWrapper(double *, unsigned long int , unsigned long int ,
					double *, unsigned long int &,
					unsigned long int &, unsigned int , SEXP *);
	void qtscoreProcessor(double *, double *, int *, int *,
					int *, int *, double *);



#ifdef __cplusplus
}
#endif

#endif
