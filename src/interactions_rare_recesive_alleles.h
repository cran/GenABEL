//=====================================================================================
//
//       Filename:  interactions_rare_recesive_alleles.h
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


#ifndef SMV_SNP_INTERACTION_RESULTS_H
#define SMV_SNP_INTERACTION_RESULTS_H


float maximumValue(float *array, unsigned number);	

double independence_test_2x2(int m[2][2], bool is_yates, int min_expected_cut_off);
enum INDEPENDENCE_TEST_NAME_ENUM {CHI2, YATES, FISHER};
double independence_test_2x2(int *x1, int *x2, int* y, int N, INDEPENDENCE_TEST_NAME_ENUM, int min_expected_cut_off);
long int factorial(int n);
double factorial(double n_);

double chi2_test(double m[2][2]);
double chi2_test_yates(double m[2][2]);
extern "C" {
double fisher_exact_test(double m[2][2]);
}

//__________________________________________________________________________________________
//the class stores interaction between two snps
class snp_snp_interaction_results
	{
	
	public:
			snp_snp_interaction_results(unsigned window_, unsigned snp_number_);
			~snp_snp_interaction_results(void);
			
			int push_chi2(float chi2, unsigned central_snp_position, unsigned window_snp_position);


			float * get_chi2_all_window(unsigned central_snp_position);
			
			float get_max_chi2(unsigned central_snp_position);
			
			float * get_maximim_chi_for_each_central_snp(void);

			unsigned get_current_window(unsigned central_snp_position);

	private:
			unsigned snp_number, window;
			float **chi2_results;
			float *maximim_chi_for_each_central_snp;
	};



#endif
