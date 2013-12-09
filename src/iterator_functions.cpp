#include <new>
#include "Rstuff.h"
#include "iterator_functions.h"

#ifdef __cplusplus
extern "C" {
#endif

// For testing purposes, tried to move some functions+wrappers to GenABEL: failed
// Product function + wrapper
double prod(double *mydata, unsigned int size)
    {
    double prodtotal = mydata[0];
    for (unsigned int i = 1; i < size; i++)
        {
        prodtotal *= mydata[i];
        }
    return prodtotal;
    }
void prodWrapper(double *indata, unsigned long int indataHeight,
                 unsigned long int indataWidth,
                 double *outdata,
                 unsigned long int &outdataNcol, unsigned long int &outdataNrow,
                 unsigned int narg, SEXP *argList)
    {
    if (indata)
        {
        outdata[0] = prod(indata, indataHeight*indataWidth);
        }
    outdataNcol = 1;
    outdataNrow = 1;
    }


// Sum function + wrapper
double sum(double *mydata, unsigned int size, bool dropNA)
    {
    double sumtotal = 0.;
    double zero = 0;
    //Rprintf("%f\n", mydata[0]);
    for (unsigned int i = 0; i < size; i++)
        {

        if (!ISNAN(mydata[i]))
            {
            sumtotal += mydata[i];
            }
        else if (!dropNA)
            {
            return(0/zero);
            }
        }
    return sumtotal;
    }

void sumWrapper(double *indata, unsigned long int indataHeight,
                unsigned long int indataWidth,
                double *outdata,
                unsigned long int &outdataNcol, unsigned long int &outdataNrow,
                unsigned int narg, SEXP *argList)
    {
    if (indata)
        {
        bool dropNA = static_cast<bool> (INTEGER(argList[0])[0]);
        outdata[0] = sum(indata, indataHeight*indataWidth, dropNA);
        }
    outdataNcol = 1;
    outdataNrow = 1;
    }

// Sum of powers function + wrapper
double sumpower(double *mydata, unsigned int size, int power)
    {
    double sumpowertotal = 0.;
    for (unsigned int i = 0; i < size; i++)
        {
        sumpowertotal += pow(mydata[i], power);
        }
    return sumpowertotal;
    }
void sumpowerWrapper(double *indata, unsigned long int indataHeight,
                     unsigned long int indataWidth,
                     double *outdata, unsigned long int &outdataNcol,
                     unsigned long int &outdataNrow, unsigned int narg, SEXP *argList)
    {
    if (indata)
        {
        int power = INTEGER(argList[0])[0];
        outdata[0] = sumpower(indata, indataHeight*indataWidth, power);
        }
    outdataNcol = 1;
    outdataNrow = 1;
    }

void qtscore_globWrapper(double *indata, unsigned long int indataHeight,
                         unsigned long int indataWidth,
                         double *outdata, unsigned long int &outdataNcol,
                         unsigned long int &outdataNrow,	unsigned int narg, SEXP *argList)
    {
    if(indata)
        {
        // Fetching data from argList
        double *pheno = REAL(argList[0]);
        int    *Type  = INTEGER(argList[1]);
        int    *Nids  = INTEGER(argList[2]);
        int    *Nstra = INTEGER(argList[3]);
        int    *stra  = INTEGER(argList[4]);
        // Calling the qtscore_glob function
        qtscoreProcessor(indata, pheno, Type, Nids, Nstra, stra, outdata);
        }
    outdataNcol = 10;
    outdataNrow = 1;
    }

void qtscoreProcessor(double *gt, double *pheno, int *Type, int *Nids, int *Nstra, int *stra,
                      double *chi2)
    {
    int nsnps = 1;;
    int nstra = (*Nstra);
    int nids =  (*Nids);
    int type = (*Type);
    int i, j, cstr, igt=0;
    int dgt;
    double Ttotg, mx;
    double Tsg0, Tsg1, Tsg2;

    double * totg = new (std::nothrow) double [nstra];
    if (totg == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * x2 = new (std::nothrow) double [nstra];
    if (x2 == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * sumx = new (std::nothrow) double [nstra];
    if (sumx == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * sg0 = new (std::nothrow) double [nstra];
    if (sg0 == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * sg1 = new (std::nothrow) double [nstra];
    if (sg1 == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * sg2 = new (std::nothrow) double [nstra];
    if (sg2 == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * xg0 = new (std::nothrow) double [nstra];
    if (xg0 == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * xg1 = new (std::nothrow) double [nstra];
    if (xg1 == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }
    double * xg2 = new (std::nothrow) double [nstra];
    if (xg2 == NULL)
        {
        Rprintf("cannot allocate RAM");
        return;
        }

    mx = -999.99;
    for (j=0; j<nstra; j++)
        {
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

    for (i=0; i<nids; i++)
        {
        if (gt[i] != 0)
            {
            cstr = stra[i];
            dgt = gt[i] - 1;
            totg[cstr]+=1.0;
            if (dgt==0)
                {
                sg0[cstr]+=1.0;
                xg0[cstr]+=pheno[i];
                }
            else if (dgt==1)
                {
                sg1[cstr]+=1.0;
                xg1[cstr]+=pheno[i];
                }
            else if (dgt==2)
                {
                sg2[cstr]+=1.0;
                xg2[cstr]+=pheno[i];
                }
            x2[cstr] += pheno[i]*pheno[i];
            sumx[cstr] += pheno[i];
            }
        }
    Ttotg=Tsg0=Tsg1=Tsg2=0.;
    for (j=0; j<nstra; j++)
        {
        Ttotg += totg[j];
        Tsg0 += sg0[j];
        Tsg1 += sg1[j];
        Tsg2 += sg2[j];
        }
    chi2[igt+6*nsnps]=Ttotg; //What is it and why igt was never initilized? (Maksim)
    if (Ttotg == 0)
        {
        chi2[igt] = -999.99;
        chi2[igt+nsnps] = -999.99;
        chi2[igt+2*nsnps] = -999.99;
        chi2[igt+3*nsnps] = -999.99;
        chi2[igt+4*nsnps] = -999.99;
        chi2[igt+5*nsnps] = -999.99;
        chi2[igt+7*nsnps] = -999.99;
        chi2[igt+8*nsnps] = -999.99;
        chi2[igt+9*nsnps] = -999.99;
        }
    else
        {
        double v, u1, u2, m0, m1, m2, v00, v02, v11, v12, v22;
        double u0;
        double bb;
        double u;

        u0 = u1 = u2 = m0 = m1 = m2 = v00 = v02 = v11 = v12 = v22 = 0.;
        for (j=0; j<nstra; j++) if (totg[j]>0)
                {
                mx = sumx[j]/totg[j];
                bb = (x2[j]/totg[j])-mx*mx;
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
        if (Tsg0>0) m0 = m0/Tsg0;
        else m0 = 1.e-16;
        if (Tsg1>0) m1 = m1/Tsg1;
        else m1 = 1.e-16;
        if (Tsg2>0) m2 = m2/Tsg2;
        else m2 = 1.e-16;


        u = u1+2.*u2;
        v = v11+4.*v12+4.*v22;
        if (v<1.e-16)
            {
            chi2[igt]=-999.99;
            chi2[igt+3*nsnps]=-999.99;
            }
        else
            {
            chi2[igt]=u*u/v;
            if (type)
                {
                double p1 = mx+u/(Tsg1+4.*Tsg2-Ttotg*((Tsg1+2.*Tsg2)/Ttotg)*((Tsg1+2.*Tsg2)/Ttotg));
                chi2[igt+3*nsnps]=(1.-mx)*p1/((1.-p1)*mx);
                }
            else
                {
                chi2[igt+3*nsnps]=u/(Tsg1+4.*Tsg2-Ttotg*((Tsg1+2.*Tsg2)/Ttotg)*((Tsg1+2.*Tsg2)/Ttotg));
                }
            }
        chi2[igt+nsnps] = -999.99;
        chi2[igt+2*nsnps] = 1.e-16;
        chi2[igt+4*nsnps] =-999.99;
        chi2[igt+5*nsnps] = -999.99;
        chi2[igt+7*nsnps] = -999.99;
        chi2[igt+8*nsnps] = -999.99;
        chi2[igt+9*nsnps] = -999.99;
        if (v00>0.)
            {
            chi2[igt+7*nsnps] = u0/sqrt(v00);
            chi2[igt+nsnps] = u0*u0/v00;
            }
        if (v22>0.)
            {
            chi2[igt+8*nsnps] = u2/sqrt(v22);
            chi2[igt+nsnps] += u2*u2/v22;
            }
        if (v00*v22>0.)
            {
            chi2[igt+9*nsnps] = v02/sqrt(v00*v22);
            chi2[igt+nsnps] += -2.*u0*u2*v02/(v00*v22);
            chi2[igt+nsnps] = chi2[igt+nsnps]/(1.-v02*v02/(v00*v22));
            }
        if (Tsg1>0)
            {
            if (type)
                {
                chi2[igt+4*nsnps]=(1.-m0)*m1/((1.-m1)*m0);
                }
            else
                {
                chi2[igt+4*nsnps]=m1-m0;
                }
            }
        if (Tsg2>0)
            {
            if (type)
                {
                chi2[igt+5*nsnps]=(1.-m0)*m2/((1.-m2)*m0);
                }
            else
                {
                chi2[igt+5*nsnps]=m2-m0;
                }
            }
        if (Tsg1>0 && Tsg2>0)
            chi2[igt+2*nsnps] = 2.;
        else if (Tsg1>0 || Tsg2>0)
            chi2[igt+2*nsnps] = 1.;
        }

    delete [] totg;
    delete [] x2;
    delete [] sumx;
    delete [] sg0;
    delete [] sg1;
    delete [] sg2;
    delete [] xg0;
    delete [] xg1;
    delete [] xg2;



    }

#ifdef __cplusplus
    }
#endif

