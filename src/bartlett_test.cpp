//#include <memory>
#include <list>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "bartlett_test.h"


using namespace std;

extern "C" {



//_____________________________________________________________
// Bartlett test for homogeneity in variances
//http://en.wikipedia.org/wiki/Bartlett%27s_test
//retrun test chi2 or -1 in case of a problem

double bartlett_test(std::list<my_small_vector> * samples)
	{
	
	unsigned sample_amount = samples->size();
	if(sample_amount <= 1) {std::cout<<"error: bartlett_test: only one sample found\n"; return -1;}

	std::list<my_small_vector>::const_iterator it;


//	double * variances = new double[sample_amount]
//	double * n = new double[sample_amount];

	double Sp2=0;
	double sum_ni_1lnSi2 = 0; //sum( (n_i-1)*ln(Si2) )
	double sum_1__n_1=0;
	double N=0;
	double var_i=0;


	for(it=samples->begin(); it!=samples->end(); ++it)
		{
		if(it->number <= 1) {std::cout<<"error: bartlett_test: one of the sample has 1 element only\n"; return -1;}
		N += it->number;
		
		var_i = var(*it);
		if(var_i > -1.E-32 && var_i < 1.E-32) {std::cout<<"error: bartlett_test: one of the sample has too small variance\n"; return -1;}

		sum_ni_1lnSi2 += (it->number - 1)*log(var_i);
		
		sum_1__n_1 += 1./(it->number-1.);

		Sp2 += (it->number - 1)*var_i;
		} 

	Sp2 = Sp2/(N-sample_amount);


	double X2 = ((N - sample_amount)*log(Sp2) - sum_ni_1lnSi2)/(1 + (sum_1__n_1-1/(N-sample_amount))/(3*(sample_amount-1)) );

	return X2;
	}





//_____________________________________________________________


double var(my_small_vector vec)
	{
	double sum = 0;
	double mean = get_mean(vec);

	if(vec.number <= 1) {std::cout<<"error: var: sample has not more than one element\n"; exit(1);}

	for(unsigned i=0; i<vec.number ; i++)
		{
		sum += (vec.vector[i]-mean)*(vec.vector[i]-mean);
		}

	return sum/(vec.number-1);
	}

//_____________________________________________________________



//_____________________________________________________________
// returun mean of the array or exit(1) in case of problem

double get_mean(my_small_vector vec)
	{
	double mean =0;
	
	if(vec.number == 0) {std::cout<<"error: get_mean: sample does not have any element\n"; exit(1);}

	for(unsigned i=0; i<vec.number ; i++)
		{
		mean += vec.vector[i];
		}

	return mean/vec.number;
	}
//_____________________________________________________________


}
