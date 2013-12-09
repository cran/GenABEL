/*
 ****************************************************
 ** GenABEL v 1.1-9 (c) 2006 Yurii Aulchenko, EMCR **
 ****************************************************
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>

int msk[4] = {192,48,12,3};
int ofs[4] = {6,4,2,0};

double SNPHWE(int, int, int);
double CalculateRS(unsigned intnAA, unsigned int nAB, unsigned int nBA, unsigned int nBB, unsigned int nDH);

void put_snps(int *b, int *Nsnps, char *a) {
	int i,j,k,idx;
	unsigned int sto;
	//  int ofs[4] = {6,4,2,0};
	int nsnps = (*Nsnps);
	int nbytes;
	if ((nsnps % 4) == 0) nbytes = nsnps/4; else nbytes = ceil(1.*nsnps/4.);
	/*
  printf("nsnps = %i\n",nsnps);
  printf("nbytes = %i\n",nbytes);
  for (i=0;i<nsnps;i++) printf ("%i ",b[i]);
  printf("\n");
	 */
	idx = 0;
	for (i=0;i<nbytes;i++) {
		sto = 0;
		for (j=0;j<4;j++) {
			k = b[idx++] << ofs[j];
			sto = sto | k;
			if (idx>=nsnps) break;
		}
		a[i] = sto;
		/*
    printf("%i %i %c\n",i,sto,a[i]);
		 */
	}
}

void decomp(char *indata, int nids, int *gt) {
	int i,j,idx;
	char str;
	//	int msk[4] = {192,48,12,3};
	//	int ofs[4] = {6,4,2,0};
	int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	idx = 0;
	for (i=0;i<nbytes;i++) {
		str = indata[i];
		for (j=0;j<4;j++) {
			gt[idx] = str & msk[j]; 
			gt[idx++] >>= ofs[j];
			//printf("%i (%i);",idx-1,gt[idx-1]);
			if (idx>=nids) {idx=0;break;}
		}
	}
	//printf("going out of decomp...\n");
}

/*
void temp(char *indata, int *Nids, int *g) {
	int i;
	int nids = (*Nids);
	printf("nids =%i\n",nids);
	decomp(indata,nids,g);
	for (i=0;i<nids;i++) printf("%i ",g[i]);
	printf("\n");
}
 */

void get_snps_many(char *a, int *Nsnps, int *Nrows, int *b) {
	int i,j,m,idx=0;
	char str;
	//  int msk[4] = {192,48,12,3};
	//  int ofs[4] = {6,4,2,0};
	int nsnps = (*Nsnps);
	int nrows = (*Nrows);
	int nbytes;
	if ((nsnps % 4) == 0) nbytes = nsnps/4; else nbytes = ceil(1.*nsnps/4.);
	/*
  printf("nsnps = %i\n",nsnps);
  printf("nbytes = %i\n",nbytes);
  printf("nrows = %i\n",nrows);
	 */
	for (m=0;m<nrows;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			/*
    printf("adr(a) = %i\n",m*nbytes + i);
			 */
			str = a[m*nbytes + i];
			/*
    Rprintf("a = %c\n",str);
    Rprintf("a = %u\n",(unsigned short int) str);
			 */
			for (j=0;j<4;j++) {
				b[m*nsnps+idx] = str & msk[j];
				b[m*nsnps+(idx++)] >>= ofs[j];
				if (idx>=nsnps) {idx=0;break;}
			}
		}
	}
	/*
  for (i=0;i<nsnps;i++) printf ("%i ",b[i]);
  printf("\n");
	 */
}

void get_snps_many_internal(char *a, int nsnps, int nrows, int *b) {
	int i,j,m,idx=0;
	char str;
	int nbytes;
	if ((nsnps % 4) == 0) nbytes = nsnps/4; else nbytes = ceil(1.*nsnps/4.);
	for (m=0;m<nrows;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = a[m*nbytes + i];
			for (j=0;j<4;j++) {
				b[m*nsnps+idx] = str & msk[j];
				b[m*nsnps+(idx++)] >>= ofs[j];
				if (idx>=nsnps) {idx=0;break;}
			}
		}
	}
}


void sset(char *indata, int *Nsnps, int *Nids, int *outlist, int *Noutlist, char *out) {
	int i,j,m,idx=0;
	char str;
	unsigned int k;
	unsigned int sto;
	//	int msk[4] = {192,48,12,3};
	//	int ofs[4] = {6,4,2,0};
	int nsnps = (*Nsnps);
	int nids = (*Nids);
	int noutlist = (*Noutlist);
	int gt[nids];
	int outgt[noutlist];
	int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	int nbyte1;
	if ((noutlist % 4) == 0) nbyte1 = noutlist/4; else nbyte1 = ceil(1.*noutlist/4.);
	//	printf("nsnps = %i\n",nsnps);
	//	printf("nbytes = %i\n",nbytes);
	//	printf("nbyte1 = %i\n",nbyte1);
	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		// subset gt
		for (i=0;i<noutlist;i++) {
			outgt[i] = gt[outlist[i]-1];
			//			printf("%i ",outgt[i]);
		} 
		//		printf("\n");
		idx = 0;
		for (i=0;i<nbyte1;i++) {
			sto = 0;
			for (j=0;j<4;j++) {
				k = outgt[idx++] << ofs[j];
				sto = sto | k; 
				if (idx>=noutlist) break;
			}
			out[m*nbyte1+i] = sto;
		}
	}
	/*
	for (i=0;i<nsnps;i++) printf ("%i ",b[i]);
	printf("\n");
	 */
}


// TO BE REMOVED LATER ON

void snp_summary(char *indata, int *Nids, int *Nsnps, double *out) {
	int i,j,m,idx;
	int nids = (*Nids);
	int nsnps = (*Nsnps);
	int gt[nids];
	char str;
	int count[3];
	//	int msk[4] = {192,48,12,3};
	//	int ofs[4] = {6,4,2,0};
	int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	double meaids,p;
	//	printf("! nsnps = %i nids = %i nbytes=%i\n",nsnps,nids,nbytes);
	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}		count[0]=count[1]=count[2]=0.;
		p = 0.;
		for (i=0;i<nids;i++)
			if (gt[i]) {
				count[gt[i]-1]++;
				p+=(gt[i]-1);
				//				printf("%i = %i  ",i,gt[i]);
			}
		meaids = 1.*(count[0]+count[1]+count[2]);
		double chi = 0.;
		double q = 2.*meaids-p;
		if (p!=0 && q!=0) {
			double div = 1./(4*meaids);
			double expec[3] = {q*q*div,2.*p*q*div,p*p*div};
			for (i=0;i<3;i++) chi+= (1.*count[i]-expec[i])*(1.*count[i]-expec[i])/expec[i];
			//			printf("%f %f %f\n",p,q,div);
			//			printf("%f %f %f\n",expec[0],expec[1],expec[2]);
		} else {chi=0.;}
		out[m]   = meaids;
		out[(nsnps)*1+m] = meaids/nids;
		out[(nsnps)*2+m] = p/(2.*meaids);
		out[(nsnps)*3+m] = count[0];
		out[(nsnps)*4+m] = count[1];
		out[(nsnps)*5+m] = count[2];
		out[(nsnps)*6+m] = chi;
		//		printf("\n%i %f \n",m,out[m*7]);
	}
}


void snp_summary_exhwe(char *indata, unsigned int *Nids, unsigned int *Nsnps, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int gt[nids];
	char str;
	unsigned int count[3];
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	double meaids,p,pmax,qmax,maf,fmax,loglik0,loglik1,chi2lrt;
	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j];
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		count[0]=count[1]=count[2]=0.;
		p = 0.;
		for (i=0;i<nids;i++)
			if (gt[i]) {
				count[gt[i]-1]++;
				p+=(gt[i]-1);
			}
		meaids = 1.*(count[0]+count[1]+count[2]);
		out[m]   = meaids;
		out[(nsnps)*1+m] = meaids/nids;
		if (meaids>0)
			out[(nsnps)*2+m] = p/(2.*meaids);
		else
			out[(nsnps)*2+m] = 0.0;
		out[(nsnps)*3+m] = count[0];
		out[(nsnps)*4+m] = count[1];
		out[(nsnps)*5+m] = count[2];
		if (meaids>0) {
			out[(nsnps)*6+m] = SNPHWE(count[1],count[0],count[2]);
			pmax = out[(nsnps)*2+m];
			qmax = 1.-pmax;
			maf = qmax; if (pmax<qmax) maf = pmax;
			if (maf>1.e-16) {
				fmax = (4.*count[0]*count[2] - 1.*count[1]*count[1])/((2.*count[0]+1.*count[1])*(2.*count[2]+1.*count[1]));
				loglik0 = 0.;
				if (count[0]) loglik0 += 2.*count[0]*log(qmax);
				if (count[1]) loglik0 += 1.*count[1]*log(2.*qmax*pmax);
				if (count[2]) loglik0 += 2.*count[2]*log(pmax);
				loglik1 = 0.;
				if (count[0]) loglik1 += 1.*count[0]*log(qmax*qmax+qmax*pmax*fmax);
				if (count[1]) loglik1 += 1.*count[1]*log(2.*qmax*pmax*(1.-fmax));
				if (count[2]) loglik1 += 1.*count[2]*log(pmax*pmax+qmax*pmax*fmax);
				chi2lrt = 2*(loglik1-loglik0);
				out[(nsnps)*7+m] = fmax;
				out[(nsnps)*8+m] = chi2lrt;
			} else {
				out[(nsnps)*7+m] = 0.;//maf;
				out[(nsnps)*8+m] = 0.;
			}
		} else {
			out[(nsnps)*6+m] = 1.0;
		}
	}
}





// END TO BE REMOVED

void redundant(char *indata, int *Nids, int *Nsnps, double *Minc, int *outlist) {
	int i,j,k,t,t1;
	int nids = (*Nids);
	int nsnps = (*Nsnps);
	double minc = (*Minc);
	int gt0[4], gt1[4];
	double maxdism = nids*(1.-minc);
	int i4 = 4, i1 = 1, s1, s2;
	int ctg[4][4];
	int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);

	for (i=0;i<nsnps;i++) outlist[i]=0;

	for (i=0;i<nsnps-1;i++) {
		if (!outlist[i])
			for (j=i+1;j<nsnps;j++) {
				if (!outlist[j]) {
					outlist[j]=i+1;
					if (minc>=1.0) {
						for (k=0;k<nbytes;k++)
							if (indata[i*nbytes+k] != indata[j*nbytes+k]) {
								outlist[j]=0;
								break;
							}
					} else {
						for (t=0;t<4;t++) for (t1=0;t1<4;t1++) ctg[t][t1]=0;
						for (k=0;k<nbytes;k++)
							if (indata[i*nbytes+k] != indata[j*nbytes+k])
							{
								get_snps_many(indata+i*nbytes+k,&i4,&i1,gt0);
								get_snps_many(indata+j*nbytes+k,&i4,&i1,gt1);
								for (t=0;t<4;t++) ctg[gt0[t]][gt1[t]]++;
							}
							else {ctg[1][1]+=4;}
						s1 = ctg[0][1]+ctg[0][2]+ctg[0][3]+
								ctg[1][0]+ctg[1][2]+ctg[1][3]+
								ctg[2][0]+ctg[2][1]+ctg[2][3]+
								ctg[3][0]+ctg[3][1]+ctg[3][2];
						s2 = ctg[0][0]+ctg[0][1]+ctg[0][3]+
								ctg[1][0]+ctg[1][2]+ctg[1][3]+
								ctg[2][1]+ctg[2][2]+ctg[2][3]+
								ctg[3][0]+ctg[3][1]+ctg[3][2];
						if (s1 > maxdism && s2 > maxdism) outlist[j]=0;
					}
				}
			}
	}
}

void fastcc_new(char *indata, int *cc, int *Nids, int *Nsnps, double *chi2) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int gt[nids];
	unsigned int count[2][4],R,N,r1,r2,n1,n2;//,S;
	double mul, a, b, c, den;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		count[0][0]=count[0][1]=count[0][2]=count[0][3]=0;
		count[1][0]=count[1][1]=count[1][2]=count[1][3]=0;
		R=N=r1=r2=n1=n2=0;
		for (i=0;i<nids;i++) count[cc[i]][gt[i]]++;
		for (j=0;j<2;j++) for (i=1;i<4;i++) {
			N+=count[j][i];
		}
		r1 = count[1][2]; r2 = count[1][3];
		n1 = r1 + count[0][2]; n2 = r2 + count[0][3];
		R = r1 + r2 + count[1][1];
		//S = count[0][1]+count[0][2]+count[0][3];
		if (N>0 && R>0 && R<N) {
			mul = (1.*N)/(1.*R*(N-R));
			a = 1.*(n1+2.*n2);
			c = 1.*(r1+2.*r2);
			b = 1.*N*c-R*a;
			den = 1.*(N*(n1+4.*n2)-a*a);
			if (den!=0.) {chi2[m] = mul*b*b/den;} else {chi2[m]=-1.;}
			chi2[m+3*nsnps] = (count[0][1])*c/((R-c)*(count[0][2]+2.*count[0][3]));
			a = 1.*(n1+n2);
			c = 1.*(r1+r2);
			b = 1.*N*c-R*a;
			den = 1.*(N*a-a*a);
			if (den!=0.) {chi2[m+nsnps] = mul*b*b/den;} else {chi2[m+nsnps]=-1.;}
			chi2[m+4*nsnps] = (count[0][1])*c/((R-c)*(count[0][2]+count[0][3]));
			a = 1.*n2;
			c = 1.*r2;
			b = 1.*N*c-R*a;
			den = 1.*(N*a-a*a);
			if (den!=0.) {chi2[m+2*nsnps] = mul*b*b/den;} else {chi2[m+2*nsnps]=-1.;}
			chi2[m+5*nsnps] = (count[0][1]+count[0][2])*c/((R-c)*(count[0][3]));
		} else {
			chi2[m] = chi2[m+nsnps] = chi2[m+2*nsnps] = chi2[m+3*nsnps] = chi2[m+4*nsnps] = chi2[m+5*nsnps] = -1.;
		}
	}
}

void fastcc(char *indata, int *cc, int *Nids, int *Nsnps, double *chi2) {
	int i,j,m,idx;
	int nids = (*Nids);
	int nsnps = (*Nsnps);
	int gt[nids];
	char str;
	int count_cas[3],count_con[3], rt[2], ct[3];
	double ecas[3], econ[3];
	double total;
	int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);

	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		count_cas[0]=count_cas[1]=count_cas[2]=count_con[0]=count_con[1]=count_con[2]=0;
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				if (cc[i]) {
					count_cas[gt[i]-1]++;
				} else {
					count_con[gt[i]-1]++;
				}
			}
		total = count_cas[0]+count_cas[1]+count_cas[2]+count_con[0]+count_con[1]+count_con[2];
		if (total==0) {
			chi2[m] = 0.;
			chi2[m+nsnps] = 0.;
			chi2[m+2*nsnps] = 0.0001;
			chi2[m+3*nsnps] = 0.;
			chi2[m+4*nsnps] = 0.;
			chi2[m+5*nsnps] = 0.;
		} else {
			/* genotypic test */
			rt[0] = count_cas[0]+count_cas[1]+count_cas[2];
			rt[1] = count_con[0]+count_con[1]+count_con[2];
			ct[0] = count_cas[0]+count_con[0];
			ct[1] = count_cas[1]+count_con[1];
			ct[2] = count_cas[2]+count_con[2];
			for (i=0;i<3;i++) ecas[i] = (1.*rt[0]*ct[i])/(1.*total);
			for (i=0;i<3;i++) econ[i] = (1.*rt[1]*ct[i])/(1.*total);
			chi2[m+nsnps] = 0.;
			for (i=0;i<3;i++) if (ecas[i]>0.) chi2[m+nsnps]+=(1.*count_cas[i]-ecas[i])*(1.*count_cas[i]-ecas[i])/ecas[i];
			for (i=0;i<3;i++) if (econ[i]>0.) chi2[m+nsnps]+=(1.*count_con[i]-econ[i])*(1.*count_con[i]-econ[i])/econ[i];
			chi2[m+2*nsnps]=2.;
			if (ct[0]<1 || ct[1]<1 || ct[2]<1) chi2[m+2*nsnps]=1.;
			if (count_cas[0]>0 && count_con[1]>0)
				chi2[m+4*nsnps] = (1.*count_cas[1]*count_con[0])/(1.*count_cas[0]*count_con[1]);
			else
				chi2[m+4*nsnps] = 10000.0;
			if (count_cas[0]>0 && count_con[2]>0)
				chi2[m+5*nsnps] = (1.*count_cas[2]*count_con[0])/(1.*count_cas[0]*count_con[2]);
			else
				chi2[m+5*nsnps] = 10000.0;

			chi2[m+6*nsnps] = rt[0] + rt[1];
			/* allelic test */
			total *= 2.;
			rt[0] *= 2.;
			rt[1] *= 2.;
			ct[0] = 2.*ct[0] + ct[1];
			ct[1] = ct[1] + 2.*ct[2];
			count_cas[0] = 2.*count_cas[0]+count_cas[1];
			count_cas[1] = count_cas[1]+2.*count_cas[2];
			count_con[0] = 2.*count_con[0]+count_con[1];
			count_con[1] = count_con[1]+2.*count_con[2];
			for (i=0;i<2;i++) ecas[i] = (1.*rt[0]*ct[i])/(1.*total);
			for (i=0;i<2;i++) econ[i] = (1.*rt[1]*ct[i])/(1.*total);
			chi2[m] = 0.;
			for (i=0;i<2;i++) if (ecas[i]>0.) chi2[m]+=(1.*count_cas[i]-ecas[i])*(1.*count_cas[i]-ecas[i])/ecas[i];
			for (i=0;i<2;i++) if (econ[i]>0.) chi2[m]+=(1.*count_con[i]-econ[i])*(1.*count_con[i]-econ[i])/econ[i];
			if (count_cas[0]>0 && count_con[1]>0)
				chi2[m+3*nsnps] = (1.*count_cas[1]*count_con[0])/(1.*count_cas[0]*count_con[1]);
			else
				chi2[m+3*nsnps] = 10000.0;
		}
	}
}

void qtscore(char *gdata, double *pheno, int *Type, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2) 
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int type = (*Type);
	int gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	double Ttotg, mx, bb, dgt, totg[nstra], x2[nstra], sumx[nstra];
	double Tsg1, Tsg2, sg1[nstra], sg2[nstra], xg1[nstra], xg2[nstra];
	double u, v, u1, u2, v11, v12, v22,det;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	//	char chgt[nbytes];

	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			x2[j] = 0.;
			sumx[j] = 0.;
			sg1[j] = 0.;
			sg2[j] = 0.;
			xg1[j] = 0.;
			xg2[j] = 0.;
		}
		for (i=0;i<nids;i++) {
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				totg[cstr]+=1.0;
				if (dgt==1) {
					sg1[cstr]+=1.0;
					xg1[cstr]+=pheno[i];
				} else if (dgt==2) {
					sg2[cstr]+=1.0;
					xg2[cstr]+=pheno[i];
				}
				x2[cstr] += pheno[i]*pheno[i];
				sumx[cstr] += pheno[i];
			}
		}
		Ttotg=Tsg1=Tsg2=0.; 
		for (j=0;j<nstra;j++) {
			Ttotg += totg[j]; 
			Tsg1 += sg1[j];
			Tsg2 += sg2[j];
		}
		chi2[igt+6*nsnps]=Ttotg;
		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u1 = u2 = v11 = v12 = v22 = 0.;
			for (j=0;j<nstra;j++) if (totg[j]>0) {
				mx = sumx[j]/totg[j];
				if (type == 0) 
					bb = (x2[j]/totg[j])-mx*mx;
				else
					bb = mx*(1.-mx);
				u1 += (xg1[j]-sg1[j]*mx);
				u2 += (xg2[j]-sg2[j]*mx);
				v11 += bb*(sg1[j]-sg1[j]*sg1[j]/totg[j]);
				v12 += bb*(0.0-sg1[j]*sg2[j]/totg[j]);
				v22 += bb*(sg2[j]-sg2[j]*sg2[j]/totg[j]);
			}
			u = u1+2.*u2;
			v = v11+4.*v12+4.*v22;
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=u*u/v;
				chi2[igt+3*nsnps]=u/(Tsg1+2.*Tsg2);
			}
			det = v11*v22 - v12*v12;
			if (det<(1.e-16)) {
				chi2[igt+nsnps]=chi2[igt];
				chi2[igt+2*nsnps] = 1.;
				chi2[igt+4*nsnps]=chi2[igt+3*nsnps];
				chi2[igt+5*nsnps]=2.*chi2[igt+3*nsnps];
			} else {
				chi2[igt+nsnps]=(u1*u1*v22+u2*u2*v11-2.0*u1*u2*v12)/det;
				chi2[igt+4*nsnps]=u1/Tsg1;
				chi2[igt+5*nsnps]=u2/Tsg2;
				if (Tsg1>0 && Tsg2>0 && Ttotg>0) 
					chi2[igt+2*nsnps] = 2.;
				else
					chi2[igt+2*nsnps] = 1.;
			}
		}
	}
}

void egscore(char *gdata, double *pheno, int *Naxes, double *axes, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2) 
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int naxes = (*Naxes);
	int gtint[nids];
	double gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	double Ttotg, mx, bb, dgt, totg[nstra], x2[nstra], sumx[nstra];
	double Tsg1, Tsg2, sg1[nstra], sg2[nstra], xg1[nstra], xg2[nstra], gamma[nstra][naxes], saxg[nstra][naxes], sa2[nstra][naxes];
	double det, u, v, u1, u2, v11, v12, v22;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	//	char chgt[nbytes];

	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gtint);
		for (i=0;i<nids;i++) {
			gt[i] = 1.*gtint[i] - 1.0;
		}
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			x2[j] = 0.;
			sumx[j] = 0.;
			sg1[j] = 0.;
			sg2[j] = 0.;
			xg1[j] = 0.;
			xg2[j] = 0.;
		}
		for (i=0;i<naxes;i++) {
			for (j=0;j<nstra;j++) {
				gamma[j][i] = 0.;
				saxg[j][i] = 0.;
				sa2[j][i] = 0.;
			}
		}
		for (j=0;j<naxes;j++) {
			for (i=0;i<nids;i++) {
				if (gtint[i] != 0) {
					cstr = stra[i];
					dgt = gt[i] - 1.0;
					saxg[cstr][j] += dgt*axes[i+j*nids];
					sa2[cstr][j] += axes[i+j*nids]*axes[i+j*nids];
					//			Rprintf("%d %d %d %e %e %e %e\n",i,j,cstr,dgt,dgt*axes[i+j*nids],saxg[cstr][j],sa2[cstr][j]);
				}
			}
		}
		for (j=0;j<naxes;j++) {
			for (i=0;i<nstra;i++) {
				gamma[i][j] = saxg[i][j]/sa2[i][j];
				//			Rprintf("%d %d %e %e %e\n",i,j,gamma[i][j],saxg[i][j],sa2[i][j]);
			}
		}
		for (i=0;i<nids;i++) {
			for (j=0;j<naxes;j++) {
				cstr = stra[i];
				gt[i] = gt[i] - gamma[cstr][j]*axes[i+j*nids];
			}
			//		Rprintf("%d %e\n",i,gt[i]);
		}
		for (j=0;j<nstra;j++) {
			for (i=0;i<nids;i++) {
				if (gtint[i] != 0) {
					cstr = stra[i];
					totg[cstr]+=1.0;
					sg1[cstr]+=gt[i];
					sg2[cstr]+=gt[i]*gt[i];
					sumx[cstr] += pheno[i];
					x2[cstr] += pheno[i]*pheno[i];
					xg1[cstr]+=gt[i]*pheno[i];
				}
			}
		}
		Ttotg=Tsg1=Tsg2=0.; 
		for (j=0;j<nstra;j++) {
			Ttotg += totg[j]; 
			Tsg1 += sg1[j];
			Tsg2 += sg2[j];
		}
		chi2[igt+6*nsnps]=Ttotg;
		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u = v = u1 = u2 = v11 = v12 = v22 = 0.;
			for (j=0;j<nstra;j++) if (totg[j]>0) {
				mx = sumx[j]/totg[j];
				bb = (x2[j]/totg[j])-mx*mx;
				u1 += (xg1[j]-sg1[j]*mx);
				u2 += (xg2[j]-sg2[j]*mx);
				v11 += bb*(sg1[j]-sg1[j]*sg1[j]/totg[j]);
				v12 += bb*(0.0-sg1[j]*sg2[j]/totg[j]);
				v22 += bb*(sg2[j]-sg2[j]*sg2[j]/totg[j]);
			}
			for (j=0;j<nstra;j++) if (totg[j]>0) {
				mx = sumx[j]/totg[j];
				bb = (x2[j]/totg[j])-mx*mx;
				u += (xg1[j]-sg1[j]*mx);
				v += bb*(sg2[j]-sg1[j]*sg1[j]/totg[j]);
			}
			//			u = u1+2.*u2;
			//			v = v11+4.*v12+4.*v22;
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=(u*u/v)*(1.*nids/(1.*nids-1.*naxes-1.));
				chi2[igt+3*nsnps]=u/(Tsg1+2.*Tsg2);
			}
			det = v11*v22 - v12*v12;
			if (det<(1.e-16)) {
				chi2[igt+nsnps]=chi2[igt];
				chi2[igt+2*nsnps] = 1.;
				chi2[igt+4*nsnps]=chi2[igt+3*nsnps];
				chi2[igt+5*nsnps]=2.*chi2[igt+3*nsnps];
			} else {
				chi2[igt+nsnps]=(u1*u1*v22+u2*u2*v11-2.0*u1*u2*v12)/det;
				chi2[igt+4*nsnps]=u1/Tsg1;
				chi2[igt+5*nsnps]=u2/Tsg2;
				if (Tsg1>0 && Tsg2>0 && Ttotg>0) 
					chi2[igt+2*nsnps] = 2.;
				else
					chi2[igt+2*nsnps] = 1.;
			}
		}
	}
}

//
//old MMSCORE
//

void mmscore(char *gdata, double *pheno, double *invS, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2) 
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	double Ttotg,dgt,totg[nstra],eG[nstra],svec[nids],gtctr[nids];
	double Tsg, sg[nstra];
	double u, v;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	//	Rprintf("nstra=%d nsnps=%d nids=%d nbytes=%d\n",nstra,nsnps,nids,nbytes);

	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		//		for (i=0;i<nids;i++) Rprintf("%d ",gt[i]-1);Rprintf("\n");
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			sg[j] = 0.;
		}
		Ttotg=Tsg=0.; 
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				totg[cstr]+=1.0;
				Ttotg += 1.0;
				sg[cstr] += dgt;
				Tsg += dgt;
			}
		chi2[igt+6*nsnps]=Ttotg;
		for (j=0;j<nstra;j++) {
			eG[j] = sg[j]/totg[j];
		}
		for (i=0;i<nids;i++) {
			gtctr[i] = 0.;
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				gtctr[i] = dgt - eG[cstr];
				//			Rprintf("i=%d gt[i]=%d gtctr[i]=%f\n",i,gt[i],gtctr[i]);
			}
		}
		for (i=0;i<nids;i++) {
			svec[i] = 0.;
			for (j=0;j<nids;j++) {
				svec[i] += gtctr[j]*invS[nids*i+j];
				//				Rprintf("%d ",nids*i+j);
			}
		}
		//		Rprintf("\nTtotg=%f\n",Ttotg);
		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u = v = 0.;
			for (i=0;i<nids;i++) {
				//				Rprintf("i=%d svec[i]=%f pheno[i]=%f gtctr[i]=%f\n",i,svec[i],pheno[i],gtctr[i]);
				if (gt[i] != 0) {
					u += svec[i]*pheno[i];
					v += svec[i]*gtctr[i];
				}
			}
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=u*u/v;
				chi2[igt+3*nsnps]=u/Tsg;
			}
			//			Rprintf("u = %f, v= %f\n",u,v);
		}
	}
}

//
//new MMSCORE (2009.01.27)
//

void mmscore_20090127(char *gdata, double *pheno, double *invS, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2) 
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	double Ttotg,dgt,totg[nstra],eG[nstra],ePH[nstra],svec[nids],gtctr[nids],phctr[nids];
	double Tsg, sg[nstra],sph[nstra];
	double u, v;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);

	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			sg[j] = 0.;
			sph[j] = 0.;
		}
		Ttotg=Tsg=0.; 
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				totg[cstr]+=1.0;
				Ttotg += 1.0;
				sg[cstr] += dgt;
				sph[cstr] += pheno[i];
				Tsg += dgt;
			}
		chi2[igt+6*nsnps]=Ttotg;
		for (j=0;j<nstra;j++) {
			eG[j] = sg[j]/totg[j];
			ePH[j] = sph[j]/totg[j];
		}
		for (i=0;i<nids;i++) {
			gtctr[i] = 0.;
			phctr[i] = 0.;
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				gtctr[i] = dgt - eG[cstr];
				phctr[i] = pheno[i] - ePH[cstr];
			}
		}
		for (i=0;i<nids;i++) {
			svec[i] = 0.;
			for (j=0;j<nids;j++) svec[i] += gtctr[j]*invS[nids*i+j];
		}
		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u = v = 0.;
			for (i=0;i<nids;i++) {
				if (gt[i] != 0) {
					u += svec[i]*phctr[i];
					v += svec[i]*gtctr[i];
				}
			}
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=u*u/v;
				chi2[igt+3*nsnps]=u/v;
			}
		}
	}
}


//
// new MMSCORE (2011.09.16)
// efficient vector-matrix multiplication
//

void mmscore_20110916(char *gdata, double *pheno, double *invS, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2)
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	double Ttotg,dgt,totg[nstra],eG[nstra],ePH[nstra],svec[nids],gtctr[nids],phctr[nids];
	double Tsg, sg[nstra],sph[nstra];
	double u, v;
	double temp1,temp2;

	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);

	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			sg[j] = 0.;
			sph[j] = 0.;
		}
		Ttotg=Tsg=0.;
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				totg[cstr]+=1.0;
				Ttotg += 1.0;
				sg[cstr] += dgt;
				sph[cstr] += pheno[i];
				Tsg += dgt;
			}
		chi2[igt+6*nsnps]=Ttotg;
		for (j=0;j<nstra;j++) {
			eG[j] = sg[j]/totg[j];
			ePH[j] = sph[j]/totg[j];
		}
		for (i=0;i<nids;i++) {
			gtctr[i] = 0.;
			phctr[i] = 0.;
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				gtctr[i] = dgt - eG[cstr];
				phctr[i] = pheno[i] - ePH[cstr];
			}
		}

		/**
		for (i=0;i<nids;i++) {
			svec[i] = 0.;
			int offS = nids*i;
			for (j=0;j<nids;j++) {
				svec[i] += gtctr[j]*invS[offS+j];
			}
		}
		 **/
		for (i=0;i<nids;i++) svec[i] = 0.;

		for (i=0;i<nids;i++) {
			int offS = nids*i;
			temp1 = gtctr[i];
			temp2 = 0.;
			for (j=(i+1);j<nids;j++) {
				svec[j] += temp1*invS[offS+j];
				temp2 += invS[offS+j]*gtctr[j];
			}
			svec[i] += temp2 + invS[offS+i]*gtctr[i];
		}

		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u = v = 0.;
			for (i=0;i<nids;i++) {
				if (gt[i] != 0) {
					u += svec[i]*phctr[i];
					v += svec[i]*gtctr[i];
				}
			}
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=u*u/v;
				chi2[igt+3*nsnps]=u/v;
			}
		}
	}
}


//
// DOES NOT WORK!!!
//

void mmscore_20110915(char *gdata, double *pheno, double *invS, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2)
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	double Ttotg,dgt,totg[nstra],eG[nstra],ePH[nstra],svec[nids],gtctr[nids];
	double Tsg, sg[nstra],sph[nstra],nIndividuals[nstra];
	double u, v;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);

	// compute strata-specific means for y
	// compute SUMi Yi SUMj OMEGAij and SUMij OMEGAij
	for (j=0;j<nstra;j++) {
		sph[j] = 0.;
		nIndividuals[j] = 0.;
		ePH[j] = 0.;
	}
	double sumOmegaJ[nids];
	double YiSumOmegaij[nids];
	double sumOmega = 0.;
//	double sumYiSumOmegaij = 0.;
	for (i=0;i<nids;i++) {
		sumOmegaJ[i] = 0.;
		YiSumOmegaij[i] = 0.;
		for (j=0;j<nids;j++) {
			cstr = stra[j];
			sph[cstr] += pheno[j];
			nIndividuals[cstr]++;
			sumOmegaJ[i] += invS[nids*i+j];
		}
		sumOmega += sumOmegaJ[i];
		YiSumOmegaij[j] = pheno[i]*sumOmegaJ[i];
//		sumYiSumOmegaij += YiSumOmegaij[j];
		ePH[j] = sph[j]/nIndividuals[j];
	}

	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			sg[j] = 0.;
		}
		Ttotg=Tsg=0.;
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				totg[cstr]+=1.0;
				Ttotg += 1.0;
				sg[cstr] += dgt;
				Tsg += dgt;
			}
		chi2[igt+6*nsnps]=Ttotg;
		for (j=0;j<nstra;j++) {
			eG[j] = sg[j]/totg[j];
		}
		for (i=0;i<nids;i++) {
			cstr = stra[i];
			if (gt[i] == 0) {
				gtctr[i] = eG[cstr];
			} else {
				gtctr[i] = 1.*gt[i] - 1.0;
			}
		}
		double sumYSumGOmega = 0.;
		double sumGSumGOmega = 0.;
//		double sumGSumOmega = 0.;
//		double SumSumGOmega = 0.;
		double summand2y, summand3y, summand2g, summand3g, summand4;
		summand2y = summand3y = summand2g = summand3g = summand4 = 0.;
		for (i=0;i<nids;i++) {
			cstr = stra[i];
			svec[i] = 0.;
			for (j=0;j<nids;j++) {
				svec[i] += gtctr[j]*invS[nids*i+j];
			}
			sumYSumGOmega += pheno[i]*svec[i];
			sumGSumGOmega += gtctr[i]*svec[i];
//			sumGSumOmega += gtctr[i]*sumOmegaJ[i];
//			SumSumGOmega += svec[i];
			summand2y += eG[cstr]*YiSumOmegaij[i];
			summand2g += eG[cstr]*gtctr[i]*sumOmegaJ[i];
			summand3y += ePH[cstr]*svec[i];
			summand3g += eG[cstr]*svec[i];
			summand4 += ePH[cstr]*eG[cstr]*sumOmega;
		}
		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u = v = 0.;
			u = sumYSumGOmega - summand2y - summand3y + summand4;
			v = sumGSumGOmega - summand2g - summand3g + summand4;
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=u*u/v;
				chi2[igt+3*nsnps]=u/v;
			}
		}
	}
}

//
// mmscore for the case of no stratification
// DOES NOT WORK (AND NO POINT)
//

void mmscore_20110915_nostrat(char *gdata, double *pheno, double *invS, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2)
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	if (nstra!=1) return;
	int nids = (*Nids);
	int gt[nids];
	int i, j, igt, i1=1;
	int nbytes;
	double Ttotg,svec[nids],gtctr[nids];//,dgt;
	double Tsg;
	double u, v;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);

	// compute strata-specific means for y
	// compute SUMi Yi SUMj OMEGAij and SUMij OMEGAij
	double OmegaX[nids*nids*4];
	double sumOmegaJ[nids];
	double YiSumOmegaij[nids];
	double sumOmega = 0.;
	double sumYiSumOmegaij = 0.;
	double sumY = 0.;
	for (i=0;i<nids;i++) {
		sumOmegaJ[i] = 0.;
		YiSumOmegaij[i] = 0.;
		for (j=0;j<nids;j++) {
			sumOmegaJ[i] += invS[nids*i+j];
			for (int omegaMult=0;omegaMult<4;omegaMult++) {
				OmegaX[nids*i+j+omegaMult*nids*nids] = invS[nids*i+j] * (1. * (omegaMult-1));
			}
		}
		sumOmega += sumOmegaJ[i];
		YiSumOmegaij[j] = pheno[i]*sumOmegaJ[i];
		sumYiSumOmegaij += YiSumOmegaij[j];
		sumY += pheno[j];
	}
	double meanY = sumY/(1.*nids);

	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		Ttotg=Tsg=0.;
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				gtctr[i] = 1.*gt[i] - 1.0;
				Ttotg += 1.0;
				Tsg += gtctr[i];
			}
		double meanG = Tsg/Ttotg;

		for (i=0;i<nids;i++) for (j=0;j<nids;j++) OmegaX[nids*i+j]=meanG;//*invS[nids*i+j];

		chi2[igt+6*nsnps]= Ttotg;
		for (i=0;i<nids;i++) {
			if (gt[i] == 0) {
				gtctr[i] = meanG;
			}
		}
		double sumYSumGOmega = 0.;
		double sumGSumGOmega = 0.;
		double sumGSumOmega = 0.;
		double SumSumGOmega = 0.;
		double summand2y, summand3y, summand2g, summand3g, summand4;
		summand2y = summand3y = summand2g = summand3g = summand4 = 0.;
//		int gtOffs[4];
//		gtOffs[0]=0;gtOffs[1]=nids*nids;gtOffs[2]=nids*nids*2;gtOffs[3]=nids*nids*3;
		for (i=0;i<nids;i++) {
			svec[i] = 0.;
			//			double * offS = &invS[nids*i];
			int offS = nids*i;
			for (j=0;j<nids;j++) {
				//				svec[i] += gtctr[j]*offS[j];
				svec[i] += gtctr[j]*invS[offS+j];
				//				svec[i] += gtctr[j]*invS[nids*i+j];
				//				svec[i] += OmegaX[offS+j+gtOffs[gt[j]]];
			}
			sumYSumGOmega += pheno[i]*svec[i];
			sumGSumGOmega += gtctr[i]*svec[i];
			sumGSumOmega += gtctr[i]*sumOmegaJ[i];
			SumSumGOmega += svec[i];
		}
		summand2g = meanG*sumGSumOmega;
		summand2y = meanG*sumYiSumOmegaij;
		summand3y = meanY*SumSumGOmega;
		summand3g = meanG*SumSumGOmega;
		summand4 = meanY*meanG*sumOmega;
		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u = v = 0.;
			u = sumYSumGOmega - summand2y - summand3y + summand4;
			v = sumGSumGOmega - summand2g - summand3g + summand4;
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=u*u/v;
				chi2[igt+3*nsnps]=u/v;
			}
		}
	}
}


void grammar(char *gdata, double *pheno, double *invS, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2) 
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	double Ttotg,dgt,totg[nstra],eG[nstra],svec[nids],gtctr[nids];
	double Tsg, sg[nstra];
	double u, v;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	//	char chgt[nbytes];

	for (i=0;i<nids;i++) {
		svec[i] = 0.;
		for (j=0;j<nids;j++) svec[i] += pheno[j]*invS[nids*j+i];
	}
	for (igt=0;igt<nsnps;igt++) {
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			sg[j] = 0.;
		}
		Ttotg=Tsg=0.; 
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				totg[cstr]+=1.0;
				Ttotg += 1.0;
				sg[cstr] += dgt;
				Tsg += dgt;
			}
		chi2[igt+6*nsnps]=Ttotg;
		for (j=0;j<nstra;j++) {
			eG[j] = sg[j]/totg[j];
		}
		for (i=0;i<nids;i++)
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = 1.*gt[i] - 1.0;
				gtctr[i] = dgt - eG[cstr];
			}
		if (Ttotg == 0) {
			chi2[igt] = 0.;
			chi2[igt+nsnps] = 0.;
			chi2[igt+2*nsnps] = 0.0001;
			chi2[igt+3*nsnps] = 0.;
			chi2[igt+4*nsnps] = 0.;
			chi2[igt+5*nsnps] = 0.;
		} else {
			u = v = 0.;
			for (i=0;i<nids;i++) {
				if (gt[i] != 0) {
					u += svec[i]*gtctr[i];
					v += gtctr[i]*gtctr[i];
				}
			}
			if (v<1.e-16) {
				chi2[igt]=0.;
				chi2[igt+3*nsnps]=0.;
			} else {
				chi2[igt]=u*u/v;
				chi2[igt+3*nsnps]=u/Tsg;
			}
		}
	}
}

void homold(char *indata, unsigned int *Nids, unsigned int *Nsnps, unsigned int *Option, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int option = (*Option);
	unsigned int gt[nids];
	unsigned int count[4],sumgt=0.;
	double homweight[4] = {0.,1.,0.,1.},p0=0.,q0=0.,maf=0.;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (i=0;i<(nids*(2+option));i++) out[i]=0.;
	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		if (option > 0) {
			count[0]=count[1]=count[2]=count[3]=sumgt=0;
			for (i=0;i<nids;i++) count[gt[i]]++;
			sumgt = count[1]+count[2]+count[3];
			p0 = (count[1]*2.+count[2]*1.)/(2.*sumgt);
			q0 = 1. - p0;
		}
		if (p0>q0) maf=q0; else maf=p0;
		for (i=0;i<nids;i++) {
			if (option == 0 && gt[i]!=0) {
				out[i]+=1.;
				out[nids+i] += homweight[gt[i]];
			} 
			if (option > 0 && gt[i]!=0 && !(maf<1.e-16) && sumgt > 1) {
				out[i]+=1.;
				out[nids+i] += homweight[gt[i]];
				out[2*nids+i] += 1. - 2.*p0*(1.-p0)*(1.*sumgt)/(1.*sumgt-1.);
			}
		}
	}
}

void hom(char *indata, unsigned int *Nids, unsigned int *Nsnps, double *freqs, double *nfreq, unsigned int *UseFreq, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int useFreq = (*UseFreq);
	unsigned int gt[nids];
	unsigned int count[4],sumgt=0.;
	double centgt[4], homweight[4] = {0.,1.,0.,1.},p0=0.,q0=0.,maf=0.;
//	double varhomweight[4] = {0.,1.,1.,1.}, centgt[4], homweight[4] = {0.,1.,0.,1.},p0=0.,q0=0.,maf=0.;
	double den, fel;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (i=0;i<(nids*5);i++) out[i]=0.;
	for (m=0;m<nsnps;m++) {
		// extract genotypes
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		// extract counts
		count[0]=count[1]=count[2]=count[3]=sumgt=0;
		for (i=0;i<nids;i++) count[gt[i]]++;
		sumgt = count[1]+count[2]+count[3];
		// extract AFs
		if (useFreq==0) {
			if (sumgt>=1) {
				p0 = (count[1]*2.+count[2]*1.)/(2.*sumgt);
				q0 = 1.-p0;
			} else {
				p0=0.;q0=1.;
			}
		} else {
			q0 = freqs[m];
			p0 = 1.-q0;
		}

		if (p0>q0) maf=q0; else maf=p0;

		if (maf>1.e-16) den=1./(p0*q0);
		else den=0;

		if (useFreq==0) {
			if (sumgt>1)
				fel =  1. - 2.*p0*q0*(1.*sumgt)/(1.*sumgt-1.);
			else
				fel = 1. - 2.*p0*q0;
		} else {
			if (nfreq[m]>1)
				fel = 1. - 2.*p0*q0*nfreq[m]/(nfreq[m]-1.);
			else
				fel = 1. - 2.*p0*q0;
		}

		centgt[0]=0.;centgt[1]=0.-q0;centgt[2]=.5-q0;centgt[3]=1.-q0;

		for (i=0;i<nids;i++)
			if (gt[i]!=0)
			{
				// compute no measured snps
				out[i]+=1.;
				// compute no measured poly snps
				if (maf>1e-16) out[nids+i]+=1.;
				// compute raw Hom
				out[2*nids+i] += homweight[gt[i]];
				// compute weighted Hom
				out[3*nids+i] += fel;
				// compute Var
				out[4*nids+i] += centgt[gt[i]]*centgt[gt[i]]*den;
				//			Rprintf("%d %d %e %e %e %d %e %e %e\n",m,i,p0,q0,maf,sumgt,nfreq[m],out[2*nids+i],1. - 2.*p0*q0*nfreq[m]/(nfreq[m]-1.));
			}

	}
}

void ibs(char *indata, unsigned int *Nids, unsigned int *Nsnps, unsigned int *Option, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int option = (*Option);
	unsigned int gt[nids],noninf=0;
	unsigned int count[4],sumgt;
	double p,q,den,centgt[4];
	double ibssum[4][4] = {{0.,0.,0.,0.},{0.,1.,.5,.0},{0.,.5,1.,.5},{0.,0.,.5,1.}};
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		noninf=0;
		if (option > 0) {
			count[0]=count[1]=count[2]=count[3]=sumgt=0;
			for (i=0;i<nids;i++) count[gt[i]]++;
			sumgt = count[1]+count[2]+count[3];
			p = (count[3]*2.+1.*count[2])/(2.*sumgt);
			q = 1.-p;
			if (p*2*sumgt < (1.-(1e-16)) || q*2*sumgt < (1.-(1e-16))) {
				noninf=1;
			} else {
				den = 1./(p*q);
				centgt[0] = 0.;
				centgt[1] = 0.-p;
				centgt[2] = .5-p;
				centgt[3] = 1.-p;
				for (i=0;i<4;i++) for (j=0;j<4;j++) ibssum[i][j] = centgt[i]*centgt[j]*den; 
			}
		}
		for (i=0;i<(nids-1);i++)
			for (j=(i+1);j<nids;j++) {
				if (gt[i]!=0 && gt[j]!=0 && !noninf) {
					out[i*nids+j]+=1.;
					out[j*nids+i]+=ibssum[gt[i]][gt[j]];
				}
			}
	}
	for (i=0;i<(nids-1);i++)
		for (j=(i+1);j<nids;j++) {
			if (out[i*nids+j]>0)
				out[j*nids+i]=out[j*nids+i]/(1.*out[i*nids+j]);
			else
				out[j*nids+i]=-1.;
		}
}


void ibsnew(char *indata, unsigned int *Nids, unsigned int *Nsnps,
		double *freqs, unsigned int *Option, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int option = (*Option);
	unsigned int gt[nids],noninf=0;
	//	unsigned int count[4],sumgt;
	double p,q,den,centgt[4];
	double ibssum[4][4] = {{0.,0.,0.,0.},{0.,1.,.5,.0},{0.,.5,1.,.5},{0.,0.,.5,1.}};
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);

	// loop over all SNPs
	for (m=0;m<nsnps;m++) {
		// extract genotypes for particular SNP in 'gt'
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}

		noninf=0;
		if (option == 1 || option == 3) { // weight = freq
			//			count[0]=count[1]=count[2]=count[3]=sumgt=0;
			//			for (i=0;i<nids;i++) count[gt[i]]++;
			//			sumgt = count[1]+count[2]+count[3];
			//			p = (count[3]*2.+1.*count[2])/(2.*sumgt);
			p = freqs[m];
			q = 1.-p;
			//			Rprintf("%d %d ",p,q);
			//			if (p*2*sumgt < (1.-(1e-16)) || q*2*sumgt < (1.-(1e-16))) {
			if (p < 1.e-16 || q < 1.e-16) {
				noninf=1;
			} else {

				centgt[0] = 0.;
				centgt[1] = 0.-p;
				centgt[2] = .5-p;
				centgt[3] = 1.-p;

				if (option == 1) {
					// use HWE assumption
					den = 1./(p*q);
				} else if (option == 3) {
					// compute empirical variance
					// note this option will give wrong results in parallel
					// implementation
					// in good way varG should be pre-computed and passed in
					// manner similar to 'p's
					double meanG = 0., nMeasured=0., ssG = 0., varG = 0., qGt;
					for (int iii = 0;iii<nids;iii++)
						if ( gt[iii] != 0) {
							qGt = centgt[ gt[iii] ];
							nMeasured += 1.0;
							meanG += qGt;
							ssG += qGt*qGt;
						}
					meanG /= nMeasured;
					varG = ssG/nMeasured - meanG*meanG;
					den = 1./(2.*varG);
				} else {
					Rprintf("Can not be!");
				}
				for (i=0;i<4;i++) for (j=0;j<4;j++) ibssum[i][j] = centgt[i]*centgt[j]*den; 
			}
		}

		for (i=0;i<(nids-1);i++)
			for (j=(i+1);j<nids;j++) {
				if (gt[i]!=0 && gt[j]!=0 && !noninf) {
					out[i*nids+j]+=1.;
					out[j*nids+i]+=ibssum[gt[i]][gt[j]];
				}
			}
	}
	// finished loop over all SNPs

	// go over all elements and divide by the number of informative markers
	for (i=0;i<(nids-1);i++)
		for (j=(i+1);j<nids;j++) {
			if (out[i*nids+j]>0)
				out[j*nids+i]=out[j*nids+i]/(1.*out[i*nids+j]);
			else
				out[j*nids+i]=-1.;
		}
}

void ibspar(char *indata, unsigned int *Nids, unsigned int *Nsnps, unsigned int *Nids1, unsigned int *ids1, unsigned int * Nids2, unsigned int *ids2, double *freqs, unsigned int *Option, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int nids1 = (*Nids1);
	unsigned int nids2 = (*Nids2);
	unsigned int option = (*Option);
	unsigned int gt[nids],noninf=0;
	double p,q,den,centgt[4];
	double ibssum[4][4] = {{0.,0.,0.,0.},{0.,1.,.5,.0},{0.,.5,1.,.5},{0.,0.,.5,1.}};
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m=0;m<nsnps;m++) {
		idx = 0;
		for (i=0;i<nbytes;i++) {
			str = indata[m*nbytes + i];
			for (j=0;j<4;j++) {
				gt[idx] = str & msk[j]; 
				gt[idx++] >>= ofs[j];
				if (idx>=nids) {idx=0;break;}
			}
		}
		noninf=0;
		if (option > 0) {
			p = freqs[m];
			q = 1.-p;
			if (p < 1.e-16 || q < 1.e-16) {
				noninf=1;
			} else {
				den = 1./(p*q);
				centgt[0] = 0.;
				centgt[1] = 0.-p;
				centgt[2] = .5-p;
				centgt[3] = 1.-p;
				for (i=0;i<4;i++) for (j=0;j<4;j++) ibssum[i][j] = centgt[i]*centgt[j]*den; 
			}
		}
		for (i=0;i<nids1;i++)
			for (j=0;j<nids2;j++) {
				if (gt[ids1[i]]!=0 && gt[ids2[j]]!=0 && !noninf) {
					out[i*nids2+j]+=ibssum[gt[ids1[i]]][gt[ids2[j]]];
					out[nids1*nids2+j*nids1+i]+=1.;
				}
			}
	}
	for (i=0;i<nids1;i++)
		for (j=0;j<nids2;j++) {
			if (out[nids1*nids2+j*nids1+i]>0)
				out[i*nids2+j]=out[i*nids2+j]/(1.*out[nids1*nids2+j*nids1+i]);
			else
				out[nids1*nids2+j*nids1+1]=-1.;
		}
}

void comp_qval(double *p, int *Length, double *out) {
	int length = (*Length);
	int i; 
	double sum,max=-1,minvec[length];//,min;
	for (i=0;i<length;i++) out[i]=0.;
	for (i=0;i<length;i++) {
		sum = i+1;
		out[i] = p[i]*(1.0*length)/sum;
		if (out[i]>max) max=out[i];
	}
	minvec[length-1]=out[length-1];
	for (i=(length-2);i>=0;i--) {
		if (out[i] < minvec[i+1]) minvec[i]=out[i]; else minvec[i]=minvec[i+1];
	}
	for (i=0;i<length;i++) {
		//min = max;
		if (out[i]<minvec[i]) out[i]=out[i]; else out[i]=minvec[i];
	}
}



void r2new(char *indata, unsigned int *Nids, unsigned int *Nsnps, unsigned int *Nsnps1, unsigned int *snps1, unsigned int *Nsnps2, unsigned int *snps2, double *out) {
	unsigned int i,j,m1,m2,idx;
	unsigned int nids = (*Nids);
//	unsigned int nsnps = (*Nsnps);
	unsigned int nsnps1 = (*Nsnps1);
	unsigned int nsnps2 = (*Nsnps2);
	unsigned int gt[2][nids],cgt[4][4],nAA,nAB,nBA,nBB,nDH;
	unsigned int csp = 0;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m1=0;m1<nsnps1;m1++)
		for (m2=0;m2<nsnps2;m2++) {
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[snps1[m1]*nbytes + i];
				for (j=0;j<4;j++) {
					gt[0][idx] = str & msk[j];
					gt[0][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[snps2[m2]*nbytes + i];
				for (j=0;j<4;j++) {
					gt[1][idx] = str & msk[j];
					gt[1][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			for (i=0;i<4;i++) for (j=0;j<4;j++) cgt[i][j]=0;
			for (i=0;i<nids;i++)
				cgt[gt[0][i]][gt[1][i]]++;
			nAA = 2*cgt[1][1] + cgt[1][2] + cgt[2][1];
			nAB = cgt[1][2] + 2*cgt[1][3] + cgt[2][3];
			nBA = cgt[2][1] + 2*cgt[3][1] + cgt[3][2];
			nBB = cgt[2][3] + cgt[3][2] + 2*cgt[3][3];
			nDH = 2*cgt[2][2];
			out[m2*nsnps1+m1] = (nAA+nAB+nBA+nBB+nDH)/2;
			if (out[m2*nsnps1+m1])
				out[nsnps1*nsnps2+m1*nsnps2+m2] = CalculateRS(nAA,nAB,nBA,nBB,nDH);
			else
				out[nsnps1*nsnps2+m1*nsnps2+m2] = 0.;
			csp++;
		}
}

void r2(char *indata, unsigned int *Nids, unsigned int *Nsnps, double *out) {
	unsigned int i,j,m0,m1,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int gt[2][nids],cgt[4][4],nAA,nAB,nBA,nBB,nDH;
	unsigned int csp = 0;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m0=0;m0<(nsnps-1);m0++)
		for (m1=(m0+1);m1<nsnps;m1++) {
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m0*nbytes + i];
				for (j=0;j<4;j++) {
					gt[0][idx] = str & msk[j];
					gt[0][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m1*nbytes + i];
				for (j=0;j<4;j++) {
					gt[1][idx] = str & msk[j];
					gt[1][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			for (i=0;i<4;i++) for (j=0;j<4;j++) cgt[i][j]=0;
			for (i=0;i<nids;i++)
				cgt[gt[0][i]][gt[1][i]]++;
			nAA = 2*cgt[1][1] + cgt[1][2] + cgt[2][1];
			nAB = cgt[1][2] + 2*cgt[1][3] + cgt[2][3];
			nBA = cgt[2][1] + 2*cgt[3][1] + cgt[3][2];
			nBB = cgt[2][3] + cgt[3][2] + 2*cgt[3][3];
			nDH = 2*cgt[2][2];
			out[m1*nsnps+m0] = (nAA+nAB+nBA+nBB+nDH)/2;
			if (out[m1*nsnps+m0])
				out[m0*nsnps+m1] = CalculateRS(nAA,nAB,nBA,nBB,nDH);
			else
				out[m0*nsnps+m1] = 0.;
			csp++;
		}
}

/**
 This code is modified code from the LdCompare program by Hao, Di and Cawley 
 (Bioinformatics, 2006, 23: 252-254)
 **/

#define R2_EM_MAX_ITERATIONS 1000
#define R2_EM_LIKELIHOOD_CONVERSION_TOLERANCE	1e-8
#define R2_EM_INITIALIZATION_FUDGE_FACTOR 0.1
#define slog(x) log((x) + 1.e-32)

double CalculateRS(unsigned int nAA, unsigned int nAB, unsigned int nBA, unsigned int nBB, unsigned int nDH)
{
	double nChrom = nAA + nAB + nBA + nBB + 2*nDH;

	/* Deal with case where either marker is monomorphic */
	short int m1_monoMorphic = (nAA+nAB == 0) || (nBA+nBB == 0);
	short int m2_monoMorphic = (nAA+nBA == 0) || (nAB+nBB == 0);
	if((m1_monoMorphic || m2_monoMorphic) && (nDH==0)) {
		return(0);
	}

	/* If we have missing data (i.e. un-phased double hets) then use EM algorithm. */
	double pAA;
	double pAB;
	double pBA;
	double pBB;
	if(nDH>0) {
		/* Set initial probs */
		double divisor = (4.0*R2_EM_INITIALIZATION_FUDGE_FACTOR) + (double) nChrom;
		pAA=(nAA + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;
		pAB=(nAB + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;
		pBA=(nBA + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;
		pBB=(nBB + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;

		double nDH_AA_BB; /* number of double hets which are AA + BB */
		double nDH_AB_BA; /* number of double hets which are AB + BA */

		double oldLogLik=-1e10;
		for(int i=0; i<R2_EM_MAX_ITERATIONS; i++) {
			/* E-step */
			double pAA_BB = pAA * pBB;
			double pAB_BA = pAB * pBA;
			nDH_AA_BB = pAA_BB/(pAA_BB + pAB_BA) * (double) nDH;
			nDH_AB_BA = ((double) nDH) - nDH_AA_BB;

			/* M-step */
			pAA = (((double)nAA) + nDH_AA_BB) / nChrom;
			pAB = (((double)nAB) + nDH_AB_BA) / nChrom;
			pBA = (((double)nBA) + nDH_AB_BA) / nChrom;
			pBB = (((double)nBB) + nDH_AA_BB) / nChrom;

			/* Iteration complete, check if we can terminate and verify that likelihood has increased */
			double logLik = ((double)nAA) * slog(pAA) + ((double)nAB) * slog(pAB) + ((double)nBA) * slog(pBA) + ((double)nBB) * slog(pBB) + ((double)nDH) * slog(pAA*pBB + pAB*pBA);
			if(i > 0) {
				//assert(logLik > oldLogLik-EPSILON); /* total likelihood should be non-decreasing, else we have a bug */
				if(logLik-oldLogLik < R2_EM_LIKELIHOOD_CONVERSION_TOLERANCE)
					break;

			}
			oldLogLik = logLik;
		}
	} else {
		pAA = ((double)nAA) / nChrom;
		pAB = ((double)nAB) / nChrom;
		pBA = ((double)nBA) / nChrom;
		pBB = ((double)nBB) / nChrom;
	}

	double p_Ax = pAA + pAB;
	double p_xA = pAA + pBA;
	double p_Bx = pBA + pBB;
	double p_xB = pAB + pBB;

	double D = pAA - p_Ax * p_xA;
	double r2 = D*D / (p_Ax * p_xA * p_Bx * p_xB);

	return r2;
}

void esthfreq(unsigned int nAA, unsigned int nAB, unsigned int nBA, unsigned int nBB, unsigned int nDH, double *eAA, double *eAB, double *eBA, double *eBB)
{
	double nChrom = nAA + nAB + nBA + nBB + 2*nDH;
	eAA[0]=eAB[0]=1;eBA[0]=eBB[0]=0;

	/* Deal with case where either marker is monomorphic */
	short int m1_monoMorphic = (nAA+nAB == 0) || (nBA+nBB == 0);
	short int m2_monoMorphic = (nAA+nBA == 0) || (nAB+nBB == 0);
	if((m1_monoMorphic || m2_monoMorphic) && (nDH==0)) {
		return;
	}

	/* If we have missing data (i.e. un-phased double hets) then use EM algorithm. */
	double pAA;
	double pAB;
	double pBA;
	double pBB;
	if(nDH>0) {
		/* Set initial probs */
		double divisor = (4.0*R2_EM_INITIALIZATION_FUDGE_FACTOR) + (double) nChrom;
		pAA=(nAA + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;
		pAB=(nAB + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;
		pBA=(nBA + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;
		pBB=(nBB + R2_EM_INITIALIZATION_FUDGE_FACTOR) / divisor;

		double nDH_AA_BB; /* number of double hets which are AA + BB */
		double nDH_AB_BA; /* number of double hets which are AB + BA */

		double oldLogLik=-1e10;
		for(int i=0; i<R2_EM_MAX_ITERATIONS; i++) {
			/* E-step */
			double pAA_BB = pAA * pBB;
			double pAB_BA = pAB * pBA;
			nDH_AA_BB = pAA_BB/(pAA_BB + pAB_BA) * (double) nDH;
			nDH_AB_BA = ((double) nDH) - nDH_AA_BB;

			/* M-step */
			pAA = (((double)nAA) + nDH_AA_BB) / nChrom;
			pAB = (((double)nAB) + nDH_AB_BA) / nChrom;
			pBA = (((double)nBA) + nDH_AB_BA) / nChrom;
			pBB = (((double)nBB) + nDH_AA_BB) / nChrom;

			/* Iteration complete, check if we can terminate and verify that likelihood has increased */
			double logLik = ((double)nAA) * slog(pAA) + ((double)nAB) * slog(pAB) + ((double)nBA) * slog(pBA) + ((double)nBB) * slog(pBB) + ((double)nDH) * slog(pAA*pBB + pAB*pBA);
			if(i > 0) {
				//assert(logLik > oldLogLik-EPSILON); /* total likelihood should be non-decreasing, else we have a bug */
				if(logLik-oldLogLik < R2_EM_LIKELIHOOD_CONVERSION_TOLERANCE)
					break;

			}
			oldLogLik = logLik;
		}
	} else {
		pAA = ((double)nAA) / nChrom;
		pAB = ((double)nAB) / nChrom;
		pBA = ((double)nBA) / nChrom;
		pBB = ((double)nBB) / nChrom;
	}

	eAA[0] = pAA*nChrom;
	eAB[0] = pAB*nChrom;
	eBA[0] = pBA*nChrom;
	eBB[0] = pBB*nChrom;

	return;
}

void rho(char *indata, unsigned int *Nids, unsigned int *Nsnps, double *out) {
	unsigned int i,j,m0,m1,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int gt[2][nids],cgt[4][4],nAA,nAB,nBA,nBB,nDH;
	double eAA,eAB,eBA,eBB,t;
	unsigned int csp = 0;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m0=0;m0<(nsnps-1);m0++)
		for (m1=(m0+1);m1<nsnps;m1++) {
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m0*nbytes + i];
				for (j=0;j<4;j++) {
					gt[0][idx] = str & msk[j];
					gt[0][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m1*nbytes + i];
				for (j=0;j<4;j++) {
					gt[1][idx] = str & msk[j];
					gt[1][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			for (i=0;i<4;i++) for (j=0;j<4;j++) cgt[i][j]=0;
			for (i=0;i<nids;i++)
				cgt[gt[0][i]][gt[1][i]]++;
			nAA = 2*cgt[1][1] + cgt[1][2] + cgt[2][1];
			nAB = cgt[1][2] + 2*cgt[1][3] + cgt[2][3];
			nBA = cgt[2][1] + 2*cgt[3][1] + cgt[3][2];
			nBB = cgt[2][3] + cgt[3][2] + 2*cgt[3][3];
			nDH = 2*cgt[2][2];
			double nchr = (double) (nAA+nAB+nBA+nBB+nDH);
			if (nchr>0) {
				esthfreq(nAA,nAB,nBA,nBB,nDH,&eAA,&eAB,&eBA,&eBB);
				if (eAA*eBB - eAB*eBA < 0) {t=eAA;eAA=eBA;eBA=t;t=eAB;eAB=eBB;eBB=t;}
				if (eAB>eBA) {t=eAA;eAA=eAB;eAB=t;t=eBA;eBA=eBB;eBB=t;}
				if (eAA*eBB - eAB*eBA < 0) {t=eAA;eAA=eBA;eBA=t;t=eAB;eAB=eBB;eBB=t;}
				if (eAB>eBA) {t=eAA;eAA=eAB;eAB=t;t=eBA;eBA=eBB;eBB=t;}
				out[m0*nsnps+m1] = (eAA*eBB-eAB*eBA)/((eAA+eAB)*(eAB+eBB));
				out[m1*nsnps+m0] = nchr*(eAA+eAB)*(eAB+eBB)/((eAA+eBA)*(eBA+eBB));
			} else {
				out[m0*nsnps+m1] = 0.;
				out[m1*nsnps+m0] = 0.;
			}
			csp++;
		}
}

void dprime(char *indata, unsigned int *Nids, unsigned int *Nsnps, double *out) {
	unsigned int i,j,m0,m1,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int gt[2][nids],cgt[4][4],nAA,nAB,nBA,nBB,nDH;
	//	double eAA,eAB,eBA,eBB,t,p_Ax,p_xA,p_Bx,p_xB,pA,pa,pB,pb,pAA,pAB,pBA,pBB,estD,estDp,Dmin,Dmax,pnAB;
	double eAA,eAB,eBA,eBB,p_Ax,p_xA,p_Bx,p_xB,pAA,pAB,pBA,pBB,estD,estDp,Dmin,Dmax;
	unsigned int csp = 0;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m0=0;m0<(nsnps-1);m0++)
		for (m1=(m0+1);m1<nsnps;m1++) {
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m0*nbytes + i];
				for (j=0;j<4;j++) {
					gt[0][idx] = str & msk[j];
					gt[0][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m1*nbytes + i];
				for (j=0;j<4;j++) {
					gt[1][idx] = str & msk[j];
					gt[1][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			for (i=0;i<4;i++) for (j=0;j<4;j++) cgt[i][j]=0;
			for (i=0;i<nids;i++)
				cgt[gt[0][i]][gt[1][i]]++;
			nAA = 2*cgt[1][1] + cgt[1][2] + cgt[2][1];
			nAB = cgt[1][2] + 2*cgt[1][3] + cgt[2][3];
			nBA = cgt[2][1] + 2*cgt[3][1] + cgt[3][2];
			nBB = cgt[2][3] + cgt[3][2] + 2*cgt[3][3];
			nDH = 2*cgt[2][2];
			double nchr = (double) (nAA+nAB+nBA+nBB+nDH);
			if (nchr>0) {
				esthfreq(nAA,nAB,nBA,nBB,nDH,&eAA,&eAB,&eBA,&eBB);
				pAA = eAA / nchr;
				pAB = eAB / nchr;
				pBA = eBA / nchr;
				pBB = eBB / nchr;
				p_Ax = pAA + pAB;
				p_xA = pAA + pBA;
				p_Bx = pBA + pBB;
				p_xB = pAB + pBB;
				/**
			pA = p_Ax; if (p_Ax < p_Bx) pA=p_Bx;
			pB = p_xA; if (p_xA < p_xB) pB=p_xB;
			pa = 1. - pA; pb = 1. - pB;
			pnAB = pAA;
			if (p_Ax>p_Bx && p_xA<p_xB) pnAB=pAB;
			if (p_Ax<p_Bx && p_xA>p_xB) pnAB=pBA;
			if (p_Ax<p_Bx && p_xA<p_xB) pnAB=pBB;
			Dmin = -1.*pA*pB; if (Dmin<(-1.*pa*pb)) Dmin=(-1.*pa*pb);
			Dmax = pA*pb; if (Dmax>pB*pa) Dmax=pB*pa;
			estD = pnAB - pA * pB;
			if (estD>0) {estDp = estD / Dmax;} else {estDp = estD / Dmin;}
				 **/
				estD = pAA*pBB - pAB*pBA;
				Dmin = p_Ax*p_xB; if (Dmin>p_Bx*p_xA) Dmin = p_Bx*p_xA;
				Dmax = -1.*p_Ax*p_xA; if (Dmax<(-1.*p_Bx*p_xB)) Dmax = -1.*p_Bx*p_xB;
				if (estD<0) estDp = estD/Dmax; else estDp = estD/Dmin;
				out[m0*nsnps+m1] = estDp;
				out[m1*nsnps+m0] = estD; //(nchr/2.);
			} else {
				out[m0*nsnps+m1] = 0.;
				out[m1*nsnps+m0] = 0.;
			}
			csp++;
		}
}

void allld(char *indata, unsigned int *Nids, unsigned int *Nsnps, double *out) {
	unsigned int i,j,m0,m1,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int gt[2][nids],cgt[4][4],nAA,nAB,nBA,nBB,nDH;
	double eAA,eAB,eBA,eBB,t;
	unsigned int csp = 0;
	char str;
	unsigned int nbytes;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	for (m0=0;m0<(nsnps-1);m0++)
		for (m1=(m0+1);m1<nsnps;m1++) {
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m0*nbytes + i];
				for (j=0;j<4;j++) {
					gt[0][idx] = str & msk[j];
					gt[0][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			idx = 0;
			for (i=0;i<nbytes;i++) {
				str = indata[m1*nbytes + i];
				for (j=0;j<4;j++) {
					gt[1][idx] = str & msk[j];
					gt[1][idx++] >>= ofs[j];
					if (idx>=nids) {idx=0;break;}
				}
			}
			for (i=0;i<4;i++) for (j=0;j<4;j++) cgt[i][j]=0;
			for (i=0;i<nids;i++)
				cgt[gt[0][i]][gt[1][i]]++;
			nAA = 2*cgt[1][1] + cgt[1][2] + cgt[2][1];
			nAB = cgt[1][2] + 2*cgt[1][3] + cgt[2][3];
			nBA = cgt[2][1] + 2*cgt[3][1] + cgt[3][2];
			nBB = cgt[2][3] + cgt[3][2] + 2*cgt[3][3];
			nDH = 2*cgt[2][2];
			double nchr = (double) (nAA+nAB+nBA+nBB+nDH);
			if (nchr>0) {
				esthfreq(nAA,nAB,nBA,nBB,nDH,&eAA,&eAB,&eBA,&eBB);
				if (eAA*eBB - eAB*eBA < 0) {t=eAA;eAA=eBA;eBA=t;t=eAB;eAB=eBB;eBB=t;}
				if (eAB>eBA) {t=eAA;eAA=eAB;eAB=t;t=eBA;eBA=eBB;eBB=t;}
				if (eAA*eBB - eAB*eBA < 0) {t=eAA;eAA=eBA;eBA=t;t=eAB;eAB=eBB;eBB=t;}
				if (eAB>eBA) {t=eAA;eAA=eAB;eAB=t;t=eBA;eBA=eBB;eBB=t;}
				out[m0*nsnps+m1] = (eAA*eBB-eAB*eBA)/((eAA+eAB)*(eAB+eBB));
				out[m1*nsnps+m0] = nchr*(eAA+eAB)*(eAB+eBB)/((eAA+eBA)*(eBA+eBB));
			} else {
				out[m0*nsnps+m1] = 0.;
				out[m1*nsnps+m0] = 0.;
			}
			csp++;
		}
}


void qtscore_glob(char *gdata, double *pheno, int *Type, int *Nids, int *Nsnps, int *Nstra, int *stra, double *chi2) 
{
	int nsnps = (*Nsnps);
	int nstra = (*Nstra);
	int nids = (*Nids);
	int type = (*Type);
	int gt[nids];
	int i, j, cstr, igt, i1=1;
	int nbytes;
	int dgt;
	double Ttotg, mx, bb, totg[nstra], x2[nstra], sumx[nstra];
	double Tsg0, Tsg1, Tsg2, sg0[nstra], sg1[nstra], sg2[nstra], xg0[nstra], xg1[nstra], xg2[nstra];
	double u, v, u0, u1, u2, m0, m1, m2, v00, v02, v11, v12, v22;//, det; 
	mx = -999.99;
	if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
	//	char chgt[nbytes];

	for (igt=0;igt<nsnps;igt++) {
		//static double det;
		get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
		for (j=0;j<nstra;j++) {
			totg[j] = 0.;
			x2[j] = 0.;
			sumx[j] = 0.;
			sg0[j] = 0.;
			sg1[j] = 0.;
			sg2[j] = 0.;
			xg0[j] = 0.;
			xg1[j] = 0.;
			xg2[j] = 0.;
		}
		for (i=0;i<nids;i++) {
			if (gt[i] != 0) {
				cstr = stra[i];
				dgt = gt[i] - 1;
				totg[cstr]+=1.0;
				if (dgt==0) {
					sg0[cstr]+=1.0;
					xg0[cstr]+=pheno[i];
				} else if (dgt==1) {
					sg1[cstr]+=1.0;
					xg1[cstr]+=pheno[i];
				} else if (dgt==2) {
					sg2[cstr]+=1.0;
					xg2[cstr]+=pheno[i];
				}
				x2[cstr] += pheno[i]*pheno[i];
				sumx[cstr] += pheno[i];
			}
		}
		Ttotg=Tsg0=Tsg1=Tsg2=0.; 
		for (j=0;j<nstra;j++) {
			Ttotg += totg[j]; 
			Tsg0 += sg0[j];
			Tsg1 += sg1[j];
			Tsg2 += sg2[j];
		}
		chi2[igt+6*nsnps]=Ttotg;
		if (Ttotg == 0) {
			chi2[igt] = -999.99;
			chi2[igt+nsnps] = -999.99;
			chi2[igt+2*nsnps] = -999.99;
			chi2[igt+3*nsnps] = -999.99;
			chi2[igt+4*nsnps] = -999.99;
			chi2[igt+5*nsnps] = -999.99;
			chi2[igt+7*nsnps] = -999.99;
			chi2[igt+8*nsnps] = -999.99;
			chi2[igt+9*nsnps] = -999.99;
		} else {
			u0 = u1 = u2 = m0 = m1 = m2 = v00 = v02 = v11 = v12 = v22 = 0.;
			for (j=0;j<nstra;j++) if (totg[j]>0) {
				mx = sumx[j]/totg[j];
				//				if (type == 0)
				bb = (x2[j]/totg[j])-mx*mx;
				//				else
				//					bb = mx*(1.-mx);
				u0 += (xg0[j]-sg0[j]*mx);
				m0 += xg0[j];
				u1 += (xg1[j]-sg1[j]*mx);
				m1 += xg1[j];
				u2 += (xg2[j]-sg2[j]*mx);
				m2 += xg2[j];
				v00 += bb*(sg0[j]-sg0[j]*sg0[j]/totg[j]);
				v11 += bb*(sg1[j]-sg1[j]*sg1[j]/totg[j]);
				v12 += bb*(0.0-sg1[j]*sg2[j]/totg[j]);
				v02 += bb*(0.0-sg0[j]*sg2[j]/totg[j]);
				v22 += bb*(sg2[j]-sg2[j]*sg2[j]/totg[j]);
			}
			if (Tsg0>0) m0 = m0/Tsg0; else m0 = 1.e-16;
			if (Tsg1>0) m1 = m1/Tsg1; else m1 = 1.e-16;
			if (Tsg2>0) m2 = m2/Tsg2; else m2 = 1.e-16;
			u = u1+2.*u2;
			v = v11+4.*v12+4.*v22;
			if (v<1.e-16) {
				chi2[igt]=-999.99;
				chi2[igt+3*nsnps]=-999.99;
			} else {
				chi2[igt]=u*u/v;
				if (type) {
					double p1 = mx+u/(Tsg1+4.*Tsg2-Ttotg*((Tsg1+2.*Tsg2)/Ttotg)*((Tsg1+2.*Tsg2)/Ttotg));
					chi2[igt+3*nsnps]=(1.-mx)*p1/((1.-p1)*mx);
				} else {
					//			  	chi2[igt+3*nsnps]=(Tsg0*(m0-mx)+Tsg1*(m1-mx)+Tsg2*(m2-mx))/Ttotg;
					chi2[igt+3*nsnps]=u/(Tsg1+4.*Tsg2-Ttotg*((Tsg1+2.*Tsg2)/Ttotg)*((Tsg1+2.*Tsg2)/Ttotg));
				}
			}
			//det = v11*v22 - v12*v12;
			//			double rho2 = v12*v12/(v11*v22);
			//			if (v00 <= 0. || v11<=0. || v22<=0. || rho2<1.e-16 || abs(det)<1.e-16) {
			chi2[igt+nsnps] = -999.99;
			chi2[igt+2*nsnps] = 1.e-16;
			chi2[igt+4*nsnps] =-999.99;
			chi2[igt+5*nsnps] = -999.99;
			chi2[igt+7*nsnps] = -999.99;
			chi2[igt+8*nsnps] = -999.99;
			chi2[igt+9*nsnps] = -999.99;
			//			} else {
			if (v00>0.) {
				chi2[igt+7*nsnps] = u0/sqrt(v00);
				chi2[igt+nsnps] = u0*u0/v00;
			}
			if (v22>0.) {
				chi2[igt+8*nsnps] = u2/sqrt(v22);
				chi2[igt+nsnps] += u2*u2/v22;
			}
			if (v00*v22>0.) {
				chi2[igt+9*nsnps] = v02/sqrt(v00*v22);
				chi2[igt+nsnps] += -2.*u0*u2*v02/(v00*v22);
				chi2[igt+nsnps] = chi2[igt+nsnps]/(1.-v02*v02/(v00*v22));
			}
			//				if (v11 > 0. && v22 > 0. && v12 > 0. && rho2<1.)
			//				if (!(v00 <= 0. || v11<=0. || v22<=0. || rho2*rho2<1.e-16 || abs(det)<1.e-16))
			//				HERE IS SOMETHING WRONG -- DO THE SAME AS IN QTSCORE CORRECTION!!!
			//				if (!(v12 <= 0. || v11<=0. || v22<=0. || (rho2*rho2-1.)<1.e-16 || abs(det)<1.e-16))
			//					chi2[igt+nsnps] = (u1*u1/v11 + u2*u2/v22 - 2.*rho2*u1*u2/v12)/(1.-rho2*rho2);
			//				else
			//					chi2[igt+nsnps] = chi2[igt];
			//(u1*u1*v22+u2*u2*v11-2.0*u1*u2*v12)/det;
			if (Tsg1>0) { if (type) {
				chi2[igt+4*nsnps]=(1.-m0)*m1/((1.-m1)*m0);
			} else {
				//				 	chi2[igt+4*nsnps]=(u1/Tsg1);
				chi2[igt+4*nsnps]=m1-m0;
			}}
			if (Tsg2>0) { if (type) {
				chi2[igt+5*nsnps]=(1.-m0)*m2/((1.-m2)*m0);
			} else {
				//			  		chi2[igt+5*nsnps]=u2/Tsg2;
				chi2[igt+5*nsnps]=m2-m0;
			}}
			if (Tsg1>0 && Tsg2>0)
				chi2[igt+2*nsnps] = 2.;
			else if (Tsg1>0 || Tsg2>0)
				chi2[igt+2*nsnps] = 1.;
			//			}
		}
	}
}
