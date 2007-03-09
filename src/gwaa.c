/*
****************************************************
** GenABEL v 1.1.3 (c) 2006 Yurii Aulchenko, EMCR **
****************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int msk[4] = {192,48,12,3};
int ofs[4] = {6,4,2,0};

double SNPHWE(int, int, int);

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
			printf("%i (%i);",idx-1,gt[idx-1]);
			if (idx>=nids) {idx=0;break;}
		}
	}
	printf("going out of decomp...\n");
}

void temp(char *indata, int *Nids, int *g) {
	int i;
	int nids = (*Nids);
	printf("nids =%i\n",nids);
	decomp(indata,nids,g);
	for (i=0;i<nids;i++) printf("%i ",g[i]);
	printf("\n");
}

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
    printf("a = %c\n",str);
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
	double meaids,p;
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
		if (meaids>0) 
			out[(nsnps)*6+m] = SNPHWE(count[1],count[0],count[2]);
		else 
			out[(nsnps)*6+m] = 1.0;
	}
}


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
	unsigned int count[2][4],R,S,N,r1,r2,n1,n2;
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
		S = count[0][1]+count[0][2]+count[0][3];
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
	double det, u, v, u1, u2, v11, v12, v22;
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

void hom(char *indata, unsigned int *Nids, unsigned int *Nsnps, unsigned int *Option, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int option = (*Option);
	unsigned int gt[nids];
	unsigned int count[4],sumgt;
	double homweight[4] = {0.,1.,0.,1.};
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
		if (option > 0) {
			count[0]=count[1]=count[2]=count[3]=sumgt=0;
			for (i=0;i<nids;i++) count[gt[i]]++;
			sumgt = count[1]+count[2]+count[3];
			homweight[1] = 1.-(count[1]*2.+count[2]*1.)/(2.*sumgt);
			homweight[3] = 1.-homweight[1];
		}
		for (i=0;i<nids;i++) {
			if (gt[i]!=0) {
				out[i]+=1.;
				out[nids+i] += homweight[gt[i]];
			}
		}
	}
}

void ibs(char *indata, unsigned int *Nids, unsigned int *Nsnps, unsigned int *Option, double *out) {
	unsigned int i,j,m,idx;
	unsigned int nids = (*Nids);
	unsigned int nsnps = (*Nsnps);
	unsigned int option = (*Option);
	unsigned int gt[nids];
	unsigned int count[4],sumgt;
	double freq[4],p,q;
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
		if (option > 0) {
			count[0]=count[1]=count[2]=count[3]=sumgt=0;
			for (i=0;i<nids;i++) count[gt[i]]++;
			sumgt = count[1]+count[2]+count[3];
			for (i=0;i<4;i++) freq[i] = (1.*count[i])/(1.*sumgt);
			p = (count[1]*2+count[2])/(2*sumgt);
			q = 1.-p;
			ibssum[1][1] = 1.-p*p;
			ibssum[2][2] = 1.-2*p*q;
			ibssum[3][3] = 1.-q*q;
			ibssum[2][1] = ibssum[0][1] = 1.-p;
			ibssum[2][3] = ibssum[2][1] = 1.-q;
		}
		for (i=0;i<(nids-1);i++)
		for (j=(i+1);j<nids;j++) {
		if (gt[i]!=0 && gt[j]!=0) {
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

void comp_qval(double *p, int *Length, double *out) {
	int length = (*Length);
	int i,j;
	double sum,max=-1,min;
	for (i=0;i<length;i++) out[i]=0.;
	for (i=0;i<length;i++) {
		sum = 0.;
		for (j=0;j<length;j++) {
			if (p[j]<=p[i]) sum += 1.;
		}
		out[i] = p[i]*(1.0*length)/sum;
		if (out[i]>max) max=out[i];
	}
	for (i=0;i<length;i++) {
		min = max;
		for (j=0;j<length;j++) {
			if (p[j]>=p[i] && out[j]<min) min = out[j];
		}
		out[i] = min;
	}
}


/*
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  

//
// Written by Jan Wigginton
*/
 
double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
   {
   if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
      {
      printf("FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a"
             " negative count", obs_hets, obs_hom1, obs_hom2);
      exit(EXIT_FAILURE);
      }

   int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
   int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

   int rare_copies = 2 * obs_homr + obs_hets;
   int genotypes   = obs_hets + obs_homc + obs_homr;

   double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
   if (het_probs == NULL)
      {
      printf("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities" );
      exit(EXIT_FAILURE);
      }

   int i;
   for (i = 0; i <= rare_copies; i++)
      het_probs[i] = 0.0;

   /* start at midpoint */
   int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

   /* check to ensure that midpoint and rare alleles have same parity */
   if ((rare_copies & 1) ^ (mid & 1))
      mid++;

   int curr_hets = mid;
   int curr_homr = (rare_copies - mid) / 2;
   int curr_homc = genotypes - curr_hets - curr_homr;

   het_probs[mid] = 1.0;
   double sum = het_probs[mid];
   for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
      {
      het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
                               / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
      sum += het_probs[curr_hets - 2];

      /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
      curr_homr++;
      curr_homc++;
      }

   curr_hets = mid;
   curr_homr = (rare_copies - mid) / 2;
   curr_homc = genotypes - curr_hets - curr_homr;
   for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
      {
      het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
                            /((curr_hets + 2.0) * (curr_hets + 1.0));
      sum += het_probs[curr_hets + 2];

      /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
      curr_homr--;
      curr_homc--;
      }

   for (i = 0; i <= rare_copies; i++)
      het_probs[i] /= sum;

   /* alternate p-value calculation for p_hi/p_lo
   double p_hi = het_probs[obs_hets];
   for (i = obs_hets + 1; i <= rare_copies; i++)
     p_hi += het_probs[i];

   double p_lo = het_probs[obs_hets];
   for (i = obs_hets - 1; i >= 0; i--)
      p_lo += het_probs[i];


   double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
   */

   double p_hwe = 0.0;
   /*  p-value calculation for p_hwe  */
   for (i = 0; i <= rare_copies; i++)
      {
      if (het_probs[i] > het_probs[obs_hets])
         continue;
      p_hwe += het_probs[i];
      }

   p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

   free(het_probs);

   return p_hwe;
   }


