//#=====================================================================================
//#
//#       Filename:  dometa.cpp
//#
//#    Description:  Function for meta analysis. 
//#
//#        Version:  1.0
//#        Created:  06-July-2009
//#       Revision:  none
//#
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================




#include<iostream>
//#include<math>

#include <Rinternals.h>

#include <R.h>  // to include Rconfig.h 

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pkg", String)
// replace pkg as appropriate 
#else
#define _(String) (String)
#endif

extern "C" {
#include "dometa.h"


void dometa_c(double *beta_set1, double *beta_set2,
			 				double *sebeta_set1, double *sebeta_set2,
						 	double *lambda_set1_, double *lambda_set2_,
			 				unsigned *num, //number of betas for current cohort
							double *mbeta,
							double *mse)
//							double *mp)
{

unsigned num_el = *num;


double wt2_set1, wt2_set2, invsumwt2;


double *lambda_set1, *lambda_set2;

if(lambda_set1_ == NULL)
	{
	lambda_set1 = new double[*num];
	for(unsigned i=0 ; i<*num ; i++) {lambda_set1[i]=1;}
	}
else
	{
	lambda_set1 = lambda_set1_;
	}

if(lambda_set2_ == NULL)
	{
	lambda_set2 = new double[*num];
	for(unsigned i=0 ; i<*num ; i++) {lambda_set2[i]=1;}
	}
else
	{
	lambda_set2 = lambda_set2_;
	}


for(unsigned i=0 ; i<num_el ; i++)
	{
	//static double se_set1; 
	//static double se_set2; 

	//se_set1 = sqrt(sebeta_set1[i]*sebeta_set1[i]*lambda_set1[i]);
	//se_set2 = sqrt(sebeta_set2[i]*sebeta_set2[i]*lambda_set2[i]);
	
	wt2_set1 = 1./(sebeta_set1[i]*sebeta_set1[i]);
	wt2_set2 = 1./(sebeta_set2[i]*sebeta_set2[i]);

	invsumwt2 = 1./(wt2_set1+wt2_set2);
	mbeta[i] = (beta_set1[i]*wt2_set1 + beta_set2[i]*wt2_set2)*invsumwt2;
	mse[i] = sqrt(invsumwt2);

//	mp[i] = (mbeta[i]*mbeta[i])/(mse[i]*mse[i]); //now it is raw mp. after we shell return mp to R we have to perform 1-pchisq(mp,1)
	}

if(lambda_set1_ == NULL)
	{
	delete[] lambda_set1;
	}

if(lambda_set2_ == NULL)
	{
	delete[] lambda_set2;
	}

}
}
