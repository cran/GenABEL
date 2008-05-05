//#=====================================================================================
//#
//#       Filename:  merge.snp.data.cpp
//#
//#    Description:  Set of functions for merging two snp.data class objects.
//#
//#        Version:  1.0
//#        Created:  18-March-2008
//#       Revision:  none
//#       
//#
//#         Author:  Maksim V. Struchalin, Yurii S. Aulchenko
//#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl, i.aoultchenko@erasmusmc.nl
//#
//#=====================================================================================

#include "gtps_container.h"
#include <cstdlib>  
#include <math.h>
#include <map>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <Rinternals.h>

#include <vector>
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


extern "C" {
typedef std::map<std::string, std::string>::const_iterator CI;



//---------------------------------------------------------
class Search
	{
	public:
	Search(unsigned *snp_intersected, unsigned *ids_intersected, unsigned num_snps_intersected, unsigned num_ids_intersected);

	unsigned what_snp_is_in_set2(unsigned id_num_in_set1);
	bool is_it_snp_in_set1(unsigned);
	unsigned what_id_is_in_set1(unsigned);


	private:
	std::map<unsigned, unsigned> array_snps_intersected_keys_from_set1;
	std::map<unsigned, unsigned> array_snps_intersected_keys_from_set2;
	std::map<unsigned, unsigned> array_ids_intersected_keys_from_set2;

	std::map<unsigned, char*> array_snps_name_intersected_keys_from_set1;

	};
//---------------------------------------------------------











//---------------------------------------------------------
Search::Search(unsigned *snp_intersected, unsigned *ids_intersected, unsigned num_snps_intersected, unsigned num_ids_intersected)
{

for(unsigned i=0 ; i<num_snps_intersected ; i++)
	{
	array_snps_intersected_keys_from_set1[snp_intersected[i]] = snp_intersected[i+num_snps_intersected];
	array_snps_intersected_keys_from_set2[snp_intersected[i+num_snps_intersected]] = snp_intersected[i];
	}


for(unsigned i=0 ; i<num_ids_intersected ; i++)
	{
	array_ids_intersected_keys_from_set2[ids_intersected[i+num_ids_intersected]] = ids_intersected[i];
	}
}
//---------------------------------------------------------



//---------------------------------------------------------
unsigned Search::what_snp_is_in_set2(unsigned snp_num_in_set1)
{
if(array_snps_intersected_keys_from_set1.find(snp_num_in_set1) == array_snps_intersected_keys_from_set1.end()) return 0;
else return array_snps_intersected_keys_from_set1[snp_num_in_set1];
}
//---------------------------------------------------------






//---------------------------------------------------------
bool Search::is_it_snp_in_set1(unsigned what_num_are_we_finding)
{
if(array_snps_intersected_keys_from_set2.find(what_num_are_we_finding) == array_snps_intersected_keys_from_set2.end()) {return false;}
else {return true;}
}
//---------------------------------------------------------






//---------------------------------------------------------
unsigned Search::what_id_is_in_set1(unsigned id_num_in_set2)
{
if(array_ids_intersected_keys_from_set2.find(id_num_in_set2) == array_ids_intersected_keys_from_set2.end()) return 0;
else return array_ids_intersected_keys_from_set2[id_num_in_set2];
}
//---------------------------------------------------------







//---------------------------------------------------------
char recoding_snp_data_under_coding_and_strand(char set1_val_, gtps_container *Set1, gtps_container *Set2, unsigned long int snp_in_set1, unsigned long int snp_in_set2, std::map<char, char*>* coding_polymorphism_map_, std::map<char, char> *alleleID_reverse_, unsigned int *snp_position_error, char* snp_set1_codding_error, char* snp_set2_codding_error,  long int *error_amount_, unsigned long int *found_error_amount_snp, unsigned long int *snp_error_counter, bool user_want_to_look_at_strand, std::map<std::string, std::string> *complemetntary)		
{
//here is coding_polymorphism_map_->first = alleleID, coding_polymorphism_map_->second="name"
char coding_set1, coding_set2; 
unsigned set1_val = unsigned(set1_val_);

//std::cout<<"recod. fun1: set1_val="<<int(set1_val)<<"\n";


coding_set1 = Set1->get_coding(snp_in_set1);
char *coding_set1_name = (*coding_polymorphism_map_)[coding_set1];

CI iter=complemetntary->find(coding_set1_name);
bool two_nucl_is_complement=false;

if(iter != complemetntary->end())
{
two_nucl_is_complement=true;
}



char coding_set1_origin=coding_set1,
		 coding_set2_origin=Set2->get_coding(snp_in_set2); 


bool made_strand_flip = false;

//recode in coding in accordance with strand
if(Set1->get_strand(snp_in_set1) != Set2->get_strand(snp_in_set2) && (user_want_to_look_at_strand || two_nucl_is_complement) )
	{
	made_strand_flip = true;
	coding_set2 = alleleID_reverse_->find(Set2->get_coding(snp_in_set2))->second;
	}
else
	{
	coding_set2 = coding_set2_origin;
	}


//std::cout<<"coding1="<<(*coding_polymorphism_map_)[coding_set1]<<", coding2="<<(*coding_polymorphism_map_)[coding_set2]<<"\n";


char *coding_set2_name;



if(coding_set1 != coding_set2)
	{
	//check wether is coding really same. For example coding_set1="AG" and coding_set2="GA" <- realy it is same polymorphism.
	coding_set2_name = (*coding_polymorphism_map_)[coding_set2];
	char coding_set2_name_new[3] = {coding_set2_name[1], coding_set2_name[0], '\0'};


	if(std::string(coding_set1_name) == std::string(coding_set2_name_new))
		 { 
				switch(set1_val)
					{
					case 1:
							set1_val=3;
							break;
					case 3:
							set1_val=1;
							break;
					default:
							break;
					}

		 }
		else
			{
			if(!made_strand_flip && Set1->get_strand(snp_in_set1) != Set2->get_strand(snp_in_set2))
		 		{
				coding_set2 = alleleID_reverse_->find(Set2->get_coding(snp_in_set2))->second;
				} //flip strand if we have not done it yet	
			
			if(coding_set1 == coding_set2) {return char(set1_val);}
			coding_set2_name = (*coding_polymorphism_map_)[coding_set2];
			char coding_set2_name_new[3] = {coding_set2_name[1], coding_set2_name[0], '\0'};

	
			if(std::string(coding_set1_name) == std::string(coding_set2_name_new))	
				{
				switch(set1_val)
					{
					case 1:
							set1_val=3;
							break;
					case 3:
							set1_val=1;
							break;
					default:
							break;
					}
				}
			else
				{
				std::ostringstream error_stream;

				if(*snp_error_counter == 0)
					{
					snp_position_error[*snp_error_counter] = snp_in_set1;
					snp_set1_codding_error[*snp_error_counter] = coding_set1_origin; 
					snp_set2_codding_error[*snp_error_counter] = coding_set2_origin;
		 			(*snp_error_counter)++;
					*found_error_amount_snp = *snp_error_counter;
					}
				else if(snp_position_error[*snp_error_counter-1] != snp_in_set1)
					{
					snp_position_error[*snp_error_counter] = snp_in_set1;
					snp_set1_codding_error[*snp_error_counter] = coding_set1_origin; 
					snp_set2_codding_error[*snp_error_counter] = coding_set2_origin;
			 		(*snp_error_counter)++;
					*found_error_amount_snp = *snp_error_counter;
					}

				if(unsigned(*snp_error_counter) >= unsigned(*error_amount_))
					{	
					std::cout<<"ID:Error: 1Too many errors while merging sets (see error table). Change error_amount value to increase error-table size.\n";
					*error_amount_=-1;
					return char();
					}

				}
	

		}
	}
//std::cout<<"recod. fun: set1_val="<<int(set1_val)<<"\n";
return char(set1_val);

}
//------------------------------------------------------------







//------------------------------------------------------------
//Function merges two data sets fast. Storages it in "const char *set_mereged". This function is performed form R code.
void fast_merge_C_(char *set1, int *num_ids1_, 
						 // Yurii: replaced with double
					 				double *num_snps1_,
			 					   char *set2, int *num_ids2_, 
						 // Yurii: replaced with double
									double *num_snps2_,
									 unsigned *num_ids_intersected, 
						 // Yurii: replaced with double
									 double *num_snps_intersected_, 
									 unsigned *snp_intersected, unsigned *ids_intersected,
									 bool *replace_na_,
									 char *strand_set1, char *strand_set2,
									 char *coding_set1, char *coding_set2,
									 char *alleleID, char **alleleID_names, unsigned *alleleID_amount,
									 char *alleleID_reverse_array,
						 // Yurii: 3 x replaced with double
									 double *error_amount_,
									 double *found_error_amount_snp_,
									 double *found_error_amount_id,
									 unsigned *id_position_error, unsigned *id_snpposition_error, char *val_set1_error, char *val_set2_error,
									 unsigned *snp_position_error, char* snp_set1_codding_error, char* snp_set2_codding_error,
									 bool *user_want_to_look_at_strand_,
									 char *set_mereged)
{
bool user_want_to_look_at_strand = *user_want_to_look_at_strand_;
bool replace_na = *replace_na_;


std::map<std::string, std::string> complemetntary;
complemetntary["AT"];
complemetntary["at"];
complemetntary["TA"];
complemetntary["ta"];
complemetntary["CG"];
complemetntary["cg"];
complemetntary["GC"];
complemetntary["gc"];









//---------------------------------
std::map<char, char*> coding_polymorphism_map;

unsigned alleleID_amount_ = *alleleID_amount;
unsigned long int error_amount = (unsigned long int) (*error_amount_);
long int error_amount_tmp = (long int) (*error_amount_); 
unsigned long int snp_error_counter=0; //this is counter for function recoding_snp_data_under_coding_and_strand 
double tmp = (*num_snps_intersected_);
unsigned long int num_snps_intersected;
unsigned long int found_error_amount_snp = (unsigned long int) (*found_error_amount_snp_);
//unsigned long int found_error_amount_id = (unsigned long int) (*found_error_amount_id_);
num_snps_intersected = (unsigned long int) tmp;


for(unsigned i=0 ; i<alleleID_amount_ ; i++)
	{
	coding_polymorphism_map[alleleID[i]] = alleleID_names[i];
	}
//---------------------------------




//---------------------------------
std::map<char, char> alleleID_reverse;
for(unsigned i=0 ; i<*alleleID_amount ; i++)
	{
	alleleID_reverse[alleleID[i]] = alleleID_reverse_array[i];
	}
//---------------------------------




Search search(snp_intersected, ids_intersected, num_snps_intersected, *num_ids_intersected);


unsigned shift_i=0;


unsigned int num_ids1  = *num_ids1_,
				 num_ids2  = *num_ids2_;
tmp = (*num_snps1_);
unsigned long int num_snps1;
num_snps1 = (unsigned long int) tmp;
tmp = (*num_snps2_);
unsigned long int num_snps2;
num_snps2 = (unsigned long int) tmp;

const unsigned long int snps_sum = num_snps1 + num_snps2;
unsigned num_ids_total_in_new_set  = num_ids1 + num_ids2 - *num_ids_intersected;  // total amount of IDs
unsigned long int num_snps_total_in_new_set = snps_sum - num_snps_intersected; // total amount of SNPs





gtps_container Set1(set1, strand_set1, coding_set1, num_ids1, num_snps1);	//creat object to facilitate working with set1
gtps_container Set2(set2, strand_set2, coding_set2, num_ids2, num_snps2);	//creat object to facilitate working with set2

gtps_container Set_merged(set_mereged, num_ids_total_in_new_set, num_snps_total_in_new_set); //Here we will storage our merged set

	

char set1_val, set2_val;
unsigned errors_counter_id=0;


//Start to merge. The mearged set will in the end 
for(unsigned long int i=1 ; i<=snps_sum ; i++)
	{
	if(i <= num_snps1)
		{

			for(unsigned long int k=1 ; k<=num_ids1 ; k++)
				{
				Set_merged.set(k, i, Set1.get(k, i));
				}//all data from current snp of set1 is copied to new set. Congratulations!!!

			
			if(unsigned long int snp_in_set2=search.what_snp_is_in_set2(i)) //snp_in_set2 is number of snp in set2 which is intersected with snp with number i 
				{
				//Variable snp_in_set2 storages snp number which corresponds with snp in set1. In other words snp with number snp_in_set2 in snp2 and current snp in set1 is similar. Not clear?
				//Now we need to copy data from set2. Without intersectioning IDs!
				unsigned long int shift_k=0;
				for(unsigned k=1 ; k<=num_ids2 ; k++)
					{
					if(unsigned id_num_in_set1 = search.what_id_is_in_set1(k)) //This means that current id and current snp belong for both set
				 		{
						set1_val = Set1.get(id_num_in_set1, i);
						set2_val = Set2.get(k, snp_in_set2);
						
//						std::cout<<"main fun.1: set1_val="<<int(set1_val)<<", set2_val="<<int(set2_val)<<"\n";

							set2_val = recoding_snp_data_under_coding_and_strand(set2_val, &Set2, &Set1, snp_in_set2, i, &coding_polymorphism_map, &alleleID_reverse , snp_position_error, snp_set1_codding_error, snp_set2_codding_error, &error_amount_tmp, &found_error_amount_snp, &snp_error_counter, user_want_to_look_at_strand, &complemetntary);	
							if(error_amount_tmp == -1)
								{
								//Function "recoding_snp_data_under_coding_and_strand" changed it to -1 because amount of errors a lot and we have to finish merging.
								return;
								}


//						std::cout<<"main fun.2: set1_val="<<int(set1_val)<<", set2_val="<<int(set2_val)<<"\n";
					
						if(set1_val != set2_val)
							{
				  // condition added by Yurii: avoid too many errors when 
				  // merging sets with largely different # of genotyped SNPs
							if (!(replace_na && (set1_val==0 || set2_val==0))) {
								id_position_error[errors_counter_id]=id_num_in_set1;
								id_snpposition_error[errors_counter_id]=i;
								val_set1_error[errors_counter_id]=set1_val;
								val_set2_error[errors_counter_id]=set2_val;
								errors_counter_id++;
							}
							if(errors_counter_id >= error_amount) 
								{
								*found_error_amount_id = (double) errors_counter_id;
								std::cout<<"ID:Error: 2Too many errors while merging sets (see error table). Change error_amount value to increase error-table size.\n";
								return;
								}	

								
								//Replace NA if it is not forbiden by user
								if(replace_na)
									{	
									if( int(set1_val) == 0 && int(set2_val) != 0 )
										{
										Set_merged.set(id_num_in_set1, i, Set2.get(k,snp_in_set2));
										}
									}
							}

						shift_k++;
					 	continue;
						}

					Set_merged.set(k - shift_k + num_ids1, i, Set2.get(k,snp_in_set2));
					}
				}


		}
	else
		{
		unsigned shift_k=0;
		
		
		if(search.is_it_snp_in_set1(i-num_snps1))
	 		{
			shift_i++; 
			continue;
			} // This snp has storaged already for both sets.
		
		for(unsigned k=1 ; k<=num_ids2 ; k++)
			{
			if(unsigned num_id_in_set1 = search.what_id_is_in_set1(k))
					{
					Set_merged.set(num_id_in_set1, i-shift_i, Set2.get(k, i-num_snps1));
					shift_k++;
					}
			else	
				{
				Set_merged.set(k-shift_k+num_ids1, i-shift_i, Set2.get(k, i-num_snps1));
				}
		
			}
			
		}


	
	}//end of for cycle




*found_error_amount_id = (double) errors_counter_id;
}
//------------------------------------------------------------








} //end of extern "C" 
