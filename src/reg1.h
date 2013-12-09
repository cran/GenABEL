#include <math.h>

class linear_reg
{
public:
	mematrix<double> beta;
	mematrix<double> sebeta;
	double sigma2;

	linear_reg(regdata rdata,int verbose)
	{
		int length_beta = (rdata.X).ncol;
//		beta = new mematrix<double> (length_beta,1);
		beta.reinit(length_beta,1);
//		sebeta = new mematrix<double> (length_beta,1);
		sebeta.reinit(length_beta,1);
		mematrix<double> tX = transpose(rdata.X);
		mematrix<double> tXX = tX*(rdata.X);

		double N = tXX.get(0,0);
		if (verbose) {Rprintf("tXX:\n");tXX.print();}
		mematrix<double> tXX_i = invert(tXX);
		if (verbose) {Rprintf("tXX-1:\n");tXX_i.print();}
		mematrix<double> tXY = tX*(rdata.Y);
		if (verbose) {Rprintf("tXY:\n");tXY.print();}
		beta = tXX_i*tXY;
		if (verbose) {Rprintf("beta:\n");(beta).print();}

		// now compute residual variance
		sigma2 = 0.;
		for (int i =0;i<(rdata.Y).nrow;i++) 
			sigma2 += ((rdata.Y).get(i,0))*((rdata.Y).get(i,0));
		for (int i=0;i<length_beta;i++) 
			sigma2 -= 2. * (beta.get(i,0)) * tXY.get(i,0);
		for (int i=0;i<(length_beta);i++) 
		for (int j=0;j<(length_beta);j++) 
			sigma2 += (beta.get(i,0)) * (beta.get(j,0)) * tXX.get(i,j);
		sigma2 /= (N - double(length_beta));

		if (verbose) {Rprintf("sigma2 = %Lf\n",sigma2);}

		for (int i=0;i<(length_beta);i++)
		{	
			double value = sqrt(sigma2*tXX_i.get(i,i));
			sebeta.put(value,i,0);
		}
		if (verbose) {Rprintf("sebeta (%d):\n",sebeta.nrow);sebeta.print();}
	}
	~linear_reg()
	{
//		delete beta;
//		delete sebeta;
	}

};

class logistic_reg
{
public:
	mematrix<double> beta;
	mematrix<double> sebeta;
	double sigma2;

	logistic_reg(regdata rdata, int verbose, int maxiter)
	{
		int length_beta = (rdata.X).ncol;
//		beta = new mematrix<double> (length_beta,1);
		beta.reinit(length_beta,1);
//		sebeta = new mematrix<double> (length_beta,1);
		sebeta.reinit(length_beta,1);
		mematrix<double> W((rdata.X).nrow,1);
		mematrix<double> z((rdata.X).nrow,1);
//		mematrix<double> eMu(rdata.X->nrow,1);
		mematrix<double> tXWX(length_beta,length_beta);
		mematrix<double> tXWX_i(length_beta,length_beta); 
		mematrix<double> tXWz(length_beta,1);

		double prev = (rdata.Y).column_mean(0);
		if (prev>=1. || prev <=0.)
		{
			//fprintf(stderr,"prevalence not within (0,1)\n");
			error("prevalence not within (0,1)");
		}
		for (int i = 0;i<length_beta;i++) beta.put(0.,i,0);
		beta.put(log(prev/(1.-prev)),0,0);
		mematrix<double> tX = transpose(rdata.X);

		//double N;

		for (int iter=0;iter<maxiter;iter++) 
		{
			mematrix<double> eMu = (rdata.X)*beta;
			for (int i=0;i<eMu.nrow;i++)
			{
				double emu = eMu.get(i,0);
				double value = emu;
				double zval = 0.;
				value = exp(value)/(1.+exp(value));
				eMu.put(value,i,0);
				W.put(value*(1.-value),i,0);
				zval = emu + (1./(value*(1.-value)))*(((rdata.Y).get(i,0))-value);
				z.put(zval,i,0);
			}
//			printMatr(eMu);
//			printMatr(z);

//			printMatr(productMatrDiag(tX,W));
			mematrix<double> tmp = productMatrDiag(tX,W);
			if (verbose) {Rprintf("tXW:\n");tmp.print();}
			mematrix<double> tXWX = tmp*(rdata.X);
			//N = tXWX.get(0,0);

			if (verbose) {Rprintf("tXWX:\n");tXWX.print();}
			tXWX_i=invert(tXWX);
			if (verbose) {Rprintf("tXWX-1:\n");tXWX_i.print();}
			mematrix<double> tmp1 = productMatrDiag(tX,W);
			mematrix<double> tXWz = tmp1*z;
			if (verbose) {Rprintf("tXWz:\n");tXWz.print();}
			beta = tXWX_i*tXWz;
			if (verbose) {Rprintf("beta:\n");beta.print();}
	
		}
		sigma2 = 0.;

		for (int i=0;i<(length_beta);i++)
		{	
			double value = sqrt(tXWX_i.get(i,i));
			sebeta.put(value,i,0);
		}
		if (verbose) {Rprintf("sebeta (%d):\n",sebeta.nrow);sebeta.print();}
	}
	~logistic_reg()
	{
//		delete beta;
//		delete sebeta;
	}
};

void coxfit2(int   *maxiter,   int   *nusedx,    int   *nvarx, 
	     double *time,      int   *status,    double *covar2, 
	     double *offset,	double *weights,   int   *strata,
	     double *means,     double *beta,      double *u, 
	     double *imat2,     double loglik[2],  int   *flag, 
	     double *work,	double *eps,       double *tol_chol,
	     double *sctest);

class coxph_reg
{
public:
	mematrix<double> beta;
	mematrix<double> sebeta;
	double sigma2;
	coxph_reg(coxph_data cdata, int maxiter, double eps, double tol_chol)
	{
		sigma2 = 0;
		beta.reinit(cdata.X.nrow,1);
		sebeta.reinit(cdata.X.nrow,1);
		mematrix<double> newoffset = cdata.offset;
		newoffset = cdata.offset - (cdata.offset).column_mean(0);

		mematrix<double> means(cdata.X.nrow,1);
		beta.reinit(cdata.X.nrow,1);
		for (int i=0;i<cdata.X.nrow;i++) beta[i]=0.;
		sebeta.reinit(cdata.X.nrow,1);
		mematrix<double> u(cdata.X.nrow,1);
		mematrix<double> imat(cdata.X.nrow,cdata.X.nrow);
		double * work = new (nothrow) double[cdata.X.ncol*2+2*(cdata.X.nrow)*(cdata.X.nrow)+3*(cdata.X.nrow)];
		if (!work) {
			//perror("can not allocate work matrix");exit(1);
			error("can not allocate work matrix");
		}
		double loglik[2];
		int flag;
		double sctest=1.0;
		
		coxfit2(&maxiter,&cdata.nids,&(cdata.X).nrow,
			cdata.stime.data,cdata.sstat.data,cdata.X.data,
			newoffset.data,cdata.weights.data,cdata.strata.data,
			means.data,beta.data,u.data,
			imat.data,loglik,&flag,
			work,&eps,&tol_chol,
			&sctest);
		delete [] work;
		for (int i=0;i<cdata.X.nrow;i++) sebeta[i]=sqrt(imat.get(i,i));
	}
	~coxph_reg()
	{
//		delete beta;
//		delete sebeta;
	}
};
