//=====================================================================================
//
//       Filename:  interactions_rare_recesive_alleles.cpp
//
//    Description:  Functions for performing chi2 test for Detecting Rare Recessive Alleles. 
//
//        Version:  1.0
//        Created:  19-May-2010
//       Revision:  none
//       
//
//         Author:  Maksim V. Struchalin
//        Company:  ErasmusMC, Epidemiology, The Netherlands.
//          Email:  m.struchalin@@erasmusmc.nl
//
//=====================================================================================



#include "gtps_container.h"
#include "interactions_rare_recesive_alleles.h"

#include <Rinternals.h>
#include <Rmath.h>



#include <Rdefines.h>
#include <vector>
#include <map>
#include <iterator>

#include <R.h>  // to include Rconfig.h 

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pkg", String)
// replace pkg as appropriate 
#else
#define _(String) (String)
#endif

//#define VARIABLE_TYPE float

extern "C" {
#include "ctest.h"
}

using namespace std;


//void
//fexact(int *nrow, int *ncol, int *table, int *ldtabl,
//				       double *expect, double *percnt, double *emin, double *prt,
//							        double *pre, /* new in C : */ int *workspace,
//											       /* new arg, was const = 30*/int *mult);






//__________________________________________________________________________________________
snp_snp_interaction_results::snp_snp_interaction_results(unsigned window_, unsigned snp_number_)
	{
	window = window_;
	snp_number = snp_number_;

	chi2_results = new float*[snp_number-1];

	unsigned window_current=window;
	for(unsigned i=0 ; i<snp_number-1 ; i++)
		{
	  if( snp_number-i+1==window_current ) window_current--;
		chi2_results[i] = new float[window_current];
		}

	maximim_chi_for_each_central_snp = new float[snp_number];

	}
//__________________________________________________________________________________________


//__________________________________________________________________________________________
snp_snp_interaction_results::~snp_snp_interaction_results(void)
	{

	for(unsigned i=0 ; i<snp_number-1 ; i++)
		{
		delete[] chi2_results[i];
		}
		
	delete[] chi2_results;
	delete[] maximim_chi_for_each_central_snp;
	}
//__________________________________________________________________________________________





//__________________________________________________________________________________________
int snp_snp_interaction_results::push_chi2(float chi2, unsigned central_snp_position, unsigned window_snp_position)
{
static unsigned int window_current;
window_current = snp_number - central_snp_position;
if(window_current > window) window_current=window;

if(central_snp_position >= snp_number) {Rprintf("snp_snp_interaction_results::push_chi2: error: central_snp_position is out of bound"); return -1;}
if(window_snp_position > window_current) {Rprintf("snp_snp_interaction_results::push_chi2: error: window is out of bound"); return -1;}	
chi2_results[central_snp_position][window_snp_position] = chi2;
//std::cout<<"snp_snp_interaction_results::push_chi2 end\n";
return 0;
}
//__________________________________________________________________________________________



//__________________________________________________________________________________________
float * snp_snp_interaction_results::get_chi2_all_window(unsigned central_snp_position)
{
if(central_snp_position > snp_number) {Rprintf("snp_snp_interaction_results::push_chi2: error: central_snp_position is out of bound"); return NULL;}
return chi2_results[central_snp_position];
}
//__________________________________________________________________________________________
			


//__________________________________________________________________________________________
float snp_snp_interaction_results::get_max_chi2(unsigned central_snp_position)
{

if(central_snp_position > snp_number) {Rprintf("snp_snp_interaction_results::push_chi2: error: central_snp_position is out of bound"); return -1;}

return maximumValue(chi2_results[central_snp_position], get_current_window(central_snp_position));
}
//__________________________________________________________________________________________



//__________________________________________________________________________________________
unsigned snp_snp_interaction_results::get_current_window(unsigned central_snp_position)
{
static unsigned int window_current;
window_current = snp_number - central_snp_position - 1;
if(window_current > window) window_current=window;

return window_current;
}
//__________________________________________________________________________________________



//__________________________________________________________________________________________
//Find maximim element
float maximumValue(float *array, unsigned number)
{
static float max;

max = NA_REAL;

//find first ellement which is not NA. If there is only NAs then return NA
//_____________________________________________
for(unsigned i = 0 ; i<number; i++)
	{
	if(ISNAN(array[i])) 
		{
		continue;
		}
	else
		{ 
		max = array[i]; break;      // start with max = first element
		}

	}

if(ISNAN(max)) return max;
//_____________________________________________


for(unsigned i = 1; i<number; i++)
	{
	if(ISNAN(array[i])) continue;
	if(array[i] > max)
	max = array[i];
	}
	
return max; // return highest value in array
}
//__________________________________________________________________________________________






//__________________________________________________________________________________________
//The function returns an array (vector) where each ellement stays for maximum value of chi2 test between a snp and the snps from left part of a window.
float * snp_snp_interaction_results::get_maximim_chi_for_each_central_snp(void)
	{
	
	for(unsigned i=0 ; i<snp_number-1 ; i++)
		{
		maximim_chi_for_each_central_snp[i] = maximumValue(get_chi2_all_window(i), get_current_window(i));
		}

	return maximim_chi_for_each_central_snp;
	}
//__________________________________________________________________________________________



//________________________________________________________________________________________
// comment:  This code is provided by "Novosibirsk programming group". Some bugs have been fixed and some correction have been introduced by me.


double independence_test_2x2(double m[2][2], INDEPENDENCE_TEST_NAME_ENUM test_name, int min_expected_cut_off)
{
static double sum;
static double expected_a, expected_b, expected_c, expected_d;


//std::cout<<" min_expected_cut_off="<<min_expected_cut_off<<"\n";

//___________________________________
// If min_expected_cut_off is specifed then do chisq test if expected value is not less then min_expected_cut_off
if(min_expected_cut_off >= 0)
	{
	sum = m[0][0] + m[0][1] + m[1][0] + m[1][1];
	expected_a = (m[0][0] + m[0][1]) * (m[0][0] + m[1][0]) / sum; //(a+b)(a+c)/n
	expected_b = (m[0][0] + m[0][1]) * (m[0][1] + m[1][1]) / sum; //(a+b)(b+d)/n
	expected_c = (m[1][0] + m[1][1]) * (m[0][0] + m[1][0]) / sum; //(c+d)(a+c)/n
	expected_d = (m[1][0] + m[1][1]) * (m[0][1] + m[1][1]) / sum; //(c+d)(b+d)/n
//	std::cout<<"expected_a="<<expected_a<<", expected_b="<<expected_b<<", expected_c="<<expected_c<<", expected_d="<<expected_d<<"\n";
	
	if(!(expected_a <= min_expected_cut_off || expected_b <= min_expected_cut_off || expected_c <= min_expected_cut_off || expected_d <= min_expected_cut_off))
		{
		test_name = CHI2;
		}
	}
//___________________________________


if(test_name == CHI2)
	{
	return chi2_test(m);
	}
else if(test_name == FISHER)
	{
	sum = m[0][0] + m[1][0] + m[0][1] + m[1][1];
	return fisher_exact_test(m); //samplesize can not exeed 100000 subjects.
	}
else if(test_name == YATES)
	{
	return chi2_test_yates(m);
	}


return -1;
}	


double chi2_test(double m[2][2])
	{
	static double sum;
	static double denom;
	static double det;

		
	sum = m[0][0] + m[1][0] + m[0][1] + m[1][1];
	denom	= (m[0][0]+m[0][1])*(m[0][0]+m[1][0])*(m[1][1]+m[1][0])*(m[1][1]+m[0][1]);
  det	= m[0][0]*m[1][1] - m[0][1]*m[1][0];
	return det*det*sum/denom;
	}


double chi2_test_yates(double m[2][2])
	{
	static double sum;
	static double denom;
	static double det;
	
	sum = m[0][0] + m[1][0] + m[0][1] + m[1][1];
	denom	= (m[0][0]+m[0][1])*(m[0][0]+m[1][0])*(m[1][1]+m[1][0])*(m[1][1]+m[0][1]);
	det	  = m[0][0]*m[1][1] - m[0][1]*m[1][0] - sum/2;
	return det*det*sum/denom;
	}

extern "C"
{
double fisher_exact_test(double m[2][2])
	{
	static int dim=2;
	
	static double expect = -1;
	static double percnt = 100;
	static double emin = 0L;
	static double prt = 1;
	static int workspace = 200000;
	static int mult = 30;
	static double pval;

	//int table_vector[4] = {int(m[0][1]+0.5), int(m[1][0]+0.5), int(m[0][1]+0.5), int(m[1][1]+0.5)};
	int table_vector[4] = {int(m[0][0]+0.5), int(m[1][0]+0.5), int(m[0][1]+0.5), int(m[1][1]+0.5)};
	
	fexact(&dim, &dim, table_vector, &dim, &expect, &percnt, &emin, &prt, &pval, &workspace, &mult);

	if(pval > 1) pval = 1;
	if(pval < 0) pval = 0;
	return qchisq(pval, 1, 0, 0);
	}
}

double independence_test_2x2(int *x1, int *x2, int* y, unsigned int N, unsigned snp1_position, unsigned snp2_position, INDEPENDENCE_TEST_NAME_ENUM test_name, int min_expected_cut_off) 
	{
	static double matrix[2][2];

	matrix[0][0]=0;
	matrix[0][1]=0;
	matrix[1][0]=0;
	matrix[1][1]=0;

	for(unsigned i=0 ; i<N ; i++)
		{
		
		if(y[i] == NA_INTEGER) {continue;}
		
		if(!(y[i] ==1 || y[i]==0)) Rf_error("Trait must posses values 0 or 1");
		

		if(x1[i]!=0 && x1[i]!=1 && x1[i]!=2 && x1[i]!=3) 
			{
			Rf_error("Snp in position %i posses unxpeted vallue for id number %i.\n", snp1_position, i);
			}
		
		if(x2[i]!=0 && x2[i]!=1 && x2[i]!=2 && x2[i]!=3) 
			{
			Rf_error("Snp in position %i posses unxpeted vallue for id number %i.\n", snp2_position, i);
			}

		if(x1[i] == 0 || x2[i] == 0) continue; // NA
		
		//if ((x1[i]==0 && x2[i]==2) || (x1[i]==2 && x2[i]==0) || (x1[i]==1 && x2[i]==1))
		if((x1[i]==1 && x2[i]==3) || (x1[i]==3 && x2[i]==1) || (x1[i]==2 && x2[i]==2))
			{
			//x3 = 1;
			matrix[1][y[i]]++;
			}
		else
			{
			//x3 = 0;
			matrix[0][y[i]]++;
			}
	}


	//In the script of Fan if at least one ellement in the table is 0 than we do not make test. For me it is not clear why.
	for(int i=0 ; i<2 ; i++)
		for(int j=0 ; j<2 ; j++)
			if(matrix[i][j]==0) return NA_REAL;
		
	return independence_test_2x2(matrix, test_name, min_expected_cut_off);	
}
//________________________________________________________________________________________







//________________________________________________________________________________________
long int factorial(int n)
	{
	if (n<=1) return(1);
	else n=n*factorial(n-1);
	return(n);
	}

double factorial(double n_)
	{
	int n = int(n_ + 0.5);
	if (n<=1) return(1);
	else n=n*factorial(n-1);
	return(n);
	}

//________________________________________________________________________________________









extern "C" {




//SEXP interaction_rare_recesive_allele_C_(char *set, int *num_ids_, int *num_snps_, int *trait_binary, int *window_)
SEXP interaction_rare_recesive_allele_C_(SEXP set_, SEXP num_ids_, SEXP num_snps_, SEXP trait_binary_, SEXP window_, SEXP return_all_result_, SEXP test_name_, SEXP min_expected_cut_off_)
{



//set - gtps slot. genabel format genotypes
//num_ids_ - id number in set
//num_snps_ - snps number in set


//double x0 = REAL(guesses)[0], 

//int num_ids = *num_ids_;
//int num_snps = *num_snps_;
//int window = *window_;



unsigned num_ids       					= unsigned(INTEGER_VALUE(num_ids_));
unsigned num_snps     					= unsigned(INTEGER_VALUE(num_snps_));
int *trait_binary								= INTEGER(trait_binary_);
unsigned window        					= unsigned(INTEGER_VALUE(window_));
int min_expected_cut_off        = INTEGER_VALUE(min_expected_cut_off_); // if observed value less then min_expected_cut_off then perform fisfer or yates
bool return_all_result 					= LOGICAL_VALUE(return_all_result_); //if true then returns vector with maximum chisq and for each snp and matrix with chi2s where each row corresponds to a snp and a column corresponds to a snp with which interaction is tested with


//_________________________________________________________________
//Which correction to use in case of small expected value
const char * test_name_char = CHARACTER_VALUE(test_name_);
INDEPENDENCE_TEST_NAME_ENUM test_name;

if(string(test_name_char) == "CHI2") test_name = CHI2;
else if(string(test_name_char) == "YATES") test_name=YATES;
else if(string(test_name_char) == "FISHER") test_name=FISHER;
else Rf_error("Test \"%s\" unknown.", test_name_char);
//_________________________________________________________________



if(min_expected_cut_off >= 0 && test_name == CHI2)
	{ 
	Rprintf("warning: Parameter min_expected_cut_off is %d and Pearson's chi-square test is chosen. Ignore min_expected_cut_off.\n", min_expected_cut_off);
	min_expected_cut_off = -1;
	}
else if(min_expected_cut_off >= 0 && test_name == YATES)
	{
	Rprintf("Running Pearson's chi-square test. Perform Yates correction in case when the expected value in contingency table below %d.\n", min_expected_cut_off );
	}
else if(min_expected_cut_off >= 0 && test_name == FISHER)
	{
	Rprintf("Running Pearson's chi-square test. Perform Fisher exact test in case when the expected value in contingency table below %d.\n", min_expected_cut_off );
	}
else if(min_expected_cut_off < 0 && test_name == YATES)
	{
	Rprintf("Running Pearson's chi-square test with yates corretion.\n");
	}
else if(min_expected_cut_off < 0 && test_name == FISHER)
	{
	Rprintf("Running Fisher exact test.\n");
	}
else if(min_expected_cut_off < 0 && test_name == CHI2)
	{
	Rprintf("Running Pearson's chi-square test.\n");
	}
	



if( num_ids > 100000 && test_name == FISHER) 
	{
	Rprintf("Number of subjects is %d that exeeds the maximum posiible value 100000. Fisher exact test can not be applied. Perform chi2 test.\n", num_ids);
	test_name=CHI2;
	}	



Rprintf("Starting analysis...\n");
//Rbyte *set = RAW(set_);
char *set = (char*)(RAW(set_));


//char* set;
//PROTECT(set_ = AS_CHARACTER(set_));
//set = R_alloc(strlen(CHAR(STRING_ELT(set_, 0))), sizeof(char));
//strcpy(set, CHAR(STRING_ELT(set_, 0)));
//UNPROTECT(1);


//std::cout<<"window="<<window<<", num_snps="<<num_snps<<", num_ids="<<num_ids<<"\n";
gtps_container Set(set, NULL, NULL, num_ids, num_snps); //creat object to facilitate working with set1


//SEXP ans;
//PROTECT(ans = allocVector(REALSXP, num_ids));
//double *rans;
//rans = REAL(ans);
//
//for(unsigned i=0 ; i<num_ids ; i++)
//	{
//	//std::cout<<"Set.get("<<i<<",1)="<<int(Set.get(i,1))<<"\n";
//	rans[i] = double(Set.get(i+1,1));
//	}
//UNPROTECT(1);
//return ans;
//Rf_error("exit");


int *snp1 = new int[num_ids];
int *snp2 = new int[num_ids];


snp_snp_interaction_results results_all(window, num_snps);

unsigned step=10000;

float chi2;
unsigned window_current=window;

//in this loop the tests between two snps within a window is done. In each iteration tests between a "central" snp and snp around (withn a window) are done.
for(unsigned snp_counter=0 ; snp_counter<num_snps-1 ; snp_counter++) //enumerate eache snp in genome
	{
//	if( num_snps-snp_counter==window_current ) window_current--; //by this we decrease window when we are close to the "left end" of the snp array.
	window_current = results_all.get_current_window(snp_counter);

	//get "center" snp
	static unsigned snp1_position;
	snp1_position=snp_counter+1;
	for(unsigned id_counter=0 ; id_counter<num_ids ; id_counter++) {snp1[id_counter] = int(Set.get(id_counter+1, snp1_position));}


	//get snps around center snp and perform chi2 test between each couple
	for(unsigned window_counter=0 ; window_counter<window_current ; window_counter++) 	
		{
	//	std::cout<<"window_current="<<window_current<<", snp_counter="<<snp_counter<<", window_counter="<<window_counter<<"\n";
		
		//get a snp
		static unsigned snp2_position;
		snp2_position=snp_counter+window_counter+2; //the first snp has position number 1
		for(unsigned id_counter=0 ; id_counter<num_ids ; id_counter++) {snp2[id_counter] = int(Set.get(id_counter+1, snp2_position));}	
	
		//perform a test
		chi2 = independence_test_2x2(snp1, snp2, trait_binary, num_ids, snp1_position, snp2_position, test_name, min_expected_cut_off);


		results_all.push_chi2(chi2, snp_counter, window_counter); 
		}

	if(snp1_position % step == 0)
		{
		Rprintf("%d SNPs done\n", snp1_position);
		if(snp1_position >= step*5) step *= 5;
		}

	}

Rprintf("All %d snps are done.\n", num_snps);

float *results = results_all.get_maximim_chi_for_each_central_snp();

SEXP results_R;





if(!return_all_result)
	{
	PROTECT(results_R = allocVector(REALSXP, num_snps-1));
	double *rresults_R = REAL(results_R);

	for(unsigned i = 0; i < num_snps-1; i++)
		{
		rresults_R[i] = results[i];
		}
	}
else
	{
	
	PROTECT(results_R = allocVector(REALSXP, num_snps-1 + (num_snps-1)*window ));
	double *rresults_R = REAL(results_R);
	
	for(unsigned i = 0; i < num_snps-1; i++)
		{
		if(ISNAN(results[i])) {rresults_R[i] = NA_REAL; continue;}
		rresults_R[i] = results[i];
		}

	unsigned current_window, counter=num_snps-1;
	float *result_window;

	for(unsigned snp=0 ; snp < num_snps-1 ; snp++)
		{

		current_window = results_all.get_current_window(snp);
		result_window = results_all.get_chi2_all_window(snp);
	
		for(unsigned window_snp=0 ; window_snp < current_window; window_snp++)
			{
		  if(ISNAN(result_window[window_snp])) {rresults_R[counter] = NA_REAL; counter++; continue;};
			rresults_R[counter] = result_window[window_snp];
			counter++;
			}
		if(current_window != window)
			{
			for(unsigned window_snp=current_window ; window_snp<window ; window_snp++)
				{
				rresults_R[counter] = NA_REAL;
				counter++;
				}
			}
			
		}
	}


UNPROTECT(1);

delete [] snp1;
delete [] snp2;

return(results_R);
} //end of interaction_rare_recesive_allele_C_



//SEXP lapply(SEXP list, SEXP expr, SEXP rho)
//		{
//		R_len_t i, n = length(list);
//		SEXP ans;
//
//		if(!isNewList(list)) error("'list' must be a list");
//		if(!isEnvironment(rho)) error("'rho' should be an environment");
//		PROTECT(ans = allocVector(VECSXP, n));
//	
//
//		for(i = 0; i < n; i++) {
//				defineVar(install("x"), VECTOR_ELT(list, i), rho);
//				SET_VECTOR_ELT(ans, i, eval(expr, rho));
//		}
//
//		
//		
//		setAttrib(ans, R_NamesSymbol, getAttrib(list, R_NamesSymbol));
//		UNPROTECT(1);
//		return(ans);
//		}

//SEXP test_eval(SEXP list, SEXP expr, SEXP rho)
//{
//		R_len_t i, n = length(list);
//		SEXP ans;
//
//		if(!isNewList(list)) error("'list' must be a list");
//		if(!isEnvironment(rho)) error("'rho' should be an environment");
//		PROTECT(ans = allocVector(VECSXP, n));
//		for(i = 0; i < n; i++) {
//				defineVar(install("x"), VECTOR_ELT(list, i), rho);
//				SET_VECTOR_ELT(ans, i, eval(expr, rho));
//		}
//		setAttrib(ans, R_NamesSymbol, getAttrib(list, R_NamesSymbol));
//		UNPROTECT(1);
//		return(ans);
//}
//
//
//SEXP fun(SEXP list)
//{
//return list;
//}

//.Call("test_eval", a, quote(.Call("fun",x)), new.env())
//.Call("test_eval", a, quote(.C("fexact",x)), new.env())





//int main(){
//		int x1[5] = {-1,1,2,2,0};
//		int x2[5] = {0,1,2,2,0};
//		int y[5] = {0,1,1,0,1};
//
//		cout << fans_test_binary(x1,x2,y,5) << endl;
//}


//________________________________________________________________________________________





//void sumpowerWrapper(double *indata, unsigned long int indataSize,
//				double *outdata, unsigned long int &outdataNcol,
//				unsigned long int &outdataNrow, unsigned int narg, double *argList) {
//		if (indata) {
//				int power = static_cast<int> (argList[0]);
//				outdata[0] = sumpower(indata, indataSize, power);
//		}
//		outdataNcol = 1;
//		outdataNrow = 1;
//}









} //end of extern "C"



