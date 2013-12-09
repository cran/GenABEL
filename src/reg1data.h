#include <string>

class regdata
{
public:
	int nids;
	int ncov;
	int noutcomes;
	int ismono;
	mematrix<double> X;
	mematrix<double> Y;

	regdata(double *OY, double *OX, int *OG, int origids, int nxcol, int nycol) 
	{
		int nmissing = 0;
		for (int i=0;i<origids;i++) if (OG[i]<0) nmissing++;
		nids = origids - nmissing;
		ncov = nxcol;
		noutcomes = nycol;
		if (nids<=0) {
			return;
		}
		X.reinit(nids,(ncov+1));
		Y.reinit(nids,noutcomes);
		int cid = 0;
		for (int j=0;j<noutcomes;j++) {
			cid = 0;
			for (int i=0;i<origids;i++) 
				if (!(OG[i]<0))	{
					Y.put(OY[j*origids+i],(cid++),j);
//					Rprintf("%d=%f\n",cid,Y[j*nids+(cid-1)]);
				}
		}
		for (int j=0;j<nxcol;j++) {
			cid = 0;
			for (int i=0;i<origids;i++) 
				if (!(OG[i]<0))	X.put(OX[j*origids+i],(cid++),j);
				
		}
		cid = 0;
		for (int i=0;i<origids;i++) 
			if (!(OG[i]<0))	X.put(double(OG[i]),(cid++),nxcol);
		ismono=1;
		for (int i=1;i<nids;i++) if (X.get(i,nxcol) != X.get(i-1,nxcol))
		{
			ismono=0;
			break;
		}
			
	}
	~regdata()
	{
//		delete X;
//		delete Y;
	}
};

// compare for sort of times
int cmpfun(const void *a, const void *b)
{
	double el1 = *(double*)a;
	double el2 = *(double*)b;
	if (el1>el2) return 1;
	if (el1<el2) return -1;
//	if (el1==el2) return 0;
	return 0;
}

class coxph_data
{
public:
	int nids;
	int ncov;
	int maxiter;
	mematrix<double> weights;
	mematrix<double> stime;
	mematrix<int>    sstat;
	mematrix<double> offset;
	mematrix<int>    strata;
	mematrix<double> X;
	mematrix<int>    order;

	coxph_data(regdata regdat) 
	{
		maxiter = 0;
		nids = regdat.nids;
		ncov = regdat.ncov;
		if (regdat.noutcomes != 2)
		{
			//fprintf(stderr,"coxph_data: number of outcomes should be 2 (now: %d)\n",regdat.noutcomes);
			error("coxph_data: number of outcomes should be 2");
		}
		X.reinit(nids,ncov);		
		stime.reinit(nids,1);
		sstat.reinit(nids,1);
		weights.reinit(nids,1);
		offset.reinit(nids,1);
		strata.reinit(nids,1);
		order.reinit(nids,1);
		for (int i=0;i<nids;i++) 
		{
			stime[i] = (regdat.Y).get(i,0);
			sstat[i] = int((regdat.Y).get(i,1));
			if (sstat[i] != 1 && sstat[i]!=0) 
			{
				//fprintf(stderr,"coxph_data: status not 0/1 (right order: id, fuptime, status ...)\n");
				error("coxph_data: status not 0/1 (right order: id, fuptime, status ...)");
			}
		}
		for (int j=0;j<ncov;j++) 
		for (int i=0;i<nids;i++) 
			X.put((regdat.X).get(i,j),i,j);

		for (int i=0;i<nids;i++) 
		{
			weights[i] = 1.0;
			offset[i] = 0.0;
			strata[i] = 0;
		}
// sort by time
		double * tmptime = new (nothrow) double [nids];
		int * passed_sorted = new (nothrow) int [nids];
		for (int i=0;i<nids;i++) {tmptime[i] = stime[i];passed_sorted[i]=0;}
		qsort(tmptime,nids,sizeof(double),cmpfun);
		for (int i=0;i<nids;i++) 
		{
			int passed = 0;
			for (int j=0;j<nids;j++)
				if (tmptime[j] == stime[i]) 
				if (!passed_sorted[j])
				{
					order[i] = j;
					passed_sorted[j] = 1;
					passed = 1;
					break;
				}
			if (passed != 1) 
			{
				//fprintf(stderr,"can not recover element %d\n",i);
				error("can not recover element");
			}
		}
		delete [] tmptime;
		delete [] passed_sorted;

		stime = reorder(stime,order);
		sstat = reorder(sstat,order);
		weights = reorder(weights,order);
		strata = reorder(strata,order);
		offset = reorder(offset,order);
		X = reorder(X,order);
		X = transpose(X);
	}
	~coxph_data()
	{
//		delete X;
//		delete sstat;
//		delete stime;
//		delete weights;
//		delete offset;
//		delete strata;
//		delete order;
	}
};


