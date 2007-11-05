#include "gtps_container.h"
#include <cstdlib>  
#include <math.h>
#include <iostream>
#include <iomanip>

#include <string>
#include <fstream>
#include <sstream>

#include <stdio.h>










//STL
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



//==================================================================
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

//	std::map<unsigned, unsigned> array_ids_intersected_back;
//	std::map<unsigned, unsigned> array_ids_intersected;


	};




extern "C" {

//Class gtps_container:

//------------------------------------------------------------------


//Class Search:


//---------------------------------------------------------
Search::Search(unsigned *snp_intersected, unsigned *ids_intersected, unsigned num_snps_intersected, unsigned num_ids_intersected)
{

//std::cout<<"num_snps_intersected="<<num_snps_intersected<<"\n";
for(unsigned i=0 ; i<num_snps_intersected ; i++)
	{
	array_snps_intersected_keys_from_set1[snp_intersected[i]] = snp_intersected[i+num_snps_intersected];
	array_snps_intersected_keys_from_set2[snp_intersected[i+num_snps_intersected]] = snp_intersected[i];
	//std::cout<<"Search::Search:   i="<<i<<", snp_intersected["<<i<<"]="	<<snp_intersected[i]<<", array_snps_intersected[snp_intersected["<<i<<"]]="<<array_snps_intersected[snp_intersected[i]]<<"\n";
	}


for(unsigned i=0 ; i<num_ids_intersected ; i++)
	{
	array_ids_intersected_keys_from_set2[ids_intersected[i+num_ids_intersected]] = ids_intersected[i];
//	array_ids_intersected[ids_intersected[i]] = ids_intersected[i+num_ids_intersected];
//	array_ids_intersected_back[ids_intersected[i+num_ids_intersected]] = ids_intersected[i];
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
//std::cout<<"Search::is_it_snp_in_set1: what_num_are_we_finding="<<what_num_are_we_finding<<"\n";
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






//==================================================================
}
















extern "C" {



const char* warning_file_name = "fast_merge_R_func_warning";



//------------------------------------------------------------
//Function merges two data sets fast. Storages it in "const char *set_mereged".
void fast_merge_C_(char *set1, int *num_ids1_, const int *num_snps1_, char *set2, const int *num_ids2_, const int *num_snps2_, unsigned *num_ids_intersected, unsigned *num_snps_intersected, unsigned *snp_intersected, unsigned *ids_intersected, unsigned *whichset, char **replace_na_, char *set_mereged)
{
//for(int i=0 ; i<10 ; i++)

enum replace_na_class {YES, NO}; 
replace_na_class replace_na;



if(std::string(*replace_na_) == "yes")
	replace_na = YES;
else
	replace_na = NO;


Search search(snp_intersected, ids_intersected, *num_snps_intersected, *num_ids_intersected);


unsigned shift_i=0;

unsigned num_ids1  = *num_ids1_,
	  		 num_snps1 = *num_snps1_,
				 num_ids2  = *num_ids2_,
	 			 num_snps2 = *num_snps2_;

const unsigned snps_sum = num_snps1 + num_snps2;
unsigned num_ids_total_in_new_set  = num_ids1 + num_ids2 - *num_ids_intersected;  // total amount of IDs
unsigned num_snps_total_in_new_set = snps_sum - *num_snps_intersected; // total amount of SNPs

//std::cout<<"num_ids1="<<num_ids1<<"\n";
//std::cout<<"num_ids2="<<num_ids2<<"\n";
//std::cout<<"num_ids_intersected="<<*num_ids_intersected<<"\n";
//std::cout<<"num_ids_total_in_new_set="<<num_ids_total_in_new_set<<"\n";

gtps_container Set1(set1, num_ids1, num_snps1);	//creat object to facilitate working with set1
gtps_container Set2(set2, num_ids2, num_snps2);	//creat object to facilitate working with set2


gtps_container Set_merged(set_mereged, num_ids_total_in_new_set, num_snps_total_in_new_set); //Here we will storage our merged set

	
std::ofstream warnig_file_stream(warning_file_name, std::ios_base::app);

//Start to merge. The mearged set will in the end 
for(unsigned i=1 ; i<=snps_sum ; i++)
	{
	if(i <= num_snps1)
		{

			

//			for(unsigned k=1 ; k<=num_ids1 ; k++)
//				{
//				std::cout<<"Set_merged.get("<<k<<", "<<i<<")="<<int(Set_merged.get(k,i))<<"\n";
//				}//all data from current snp of set1 is copied to new set. Congratulations!!!




			//----------------------
			//----------------------
			for(unsigned k=1 ; k<=num_ids1 ; k++)
				{
				Set_merged.set(k, i, Set1.get(k, i));
				//std::cout<<"Set_merged.set("<<k<<", "<<i<<", "<<int(Set1.get(k, i))<<")\n";
			//	std::cout<<"Set_merged.get("<<k<<", "<<i<<")="<<int(Set_merged.get(k,i))<<"\n";
//					for(unsigned k=1 ; k<=num_ids1 ; k++)
//						{
//						std::cout<<"Set_merged.get("<<k<<", "<<i<<")="<<int(Set_merged.get(k,i))<<"\n";
//						}//all data from current snp of set1 is copied to new set. Congratulations!!!
				}//all data from current snp of set1 is copied to new set. Congratulations!!!

//if(i == 2) return ;
			

//return;

//			std::cout<<"search.what_snp_is_in_set2("<<i<<")="<<search.what_snp_is_in_set2(i)<<"\n";
			if(unsigned snp_in_set2=search.what_snp_is_in_set2(i)) //snp_in_set2 is number of snp in set2 which is intersected with snp with number i 
				{
				//Variable snp_in_set2 storages snp number which corresponds with snp in set1. In other words snp with number snp_in_set2 in snp2 and current snp in set1 is similar. Not clear?
			
				//Now we need to copy data from set2. Without intersectioning IDs! Difficult task?
				unsigned int shift_k=0;
				
				for(unsigned k=1 ; k<=num_ids2 ; k++)
					{
					if(unsigned id_num_in_set1 = search.what_id_is_in_set1(k)) //This means that current id and current snp belong for both set
				 		{
						//Now we need to check wether polymorfism in set1 and set2 (taking into account strands)
						//for current snp and id is same.


						if(Set1.get(id_num_in_set1, i) != Set2.get(k, snp_in_set2))
							{
							//std::cout<<"SNP for for id with number "<<id_num_in_set1<<" for snp with number "<<i<<" in set1 is anothers than in set2\n";
							warnig_file_stream<<"SNP for for id with number "<<id_num_in_set1<<"for snp with number "<<i<<" in set1 is anothers than in set2\n";



							if(*whichset == 1)
								{
								//If in set1 this snp for current id is not measured but in set2 is measured than replace it by data from ste2
								if( int(Set1.get(id_num_in_set1, i)) == 0 && int(Set2.get(k, snp_in_set2)) != 0 )
									{
									warnig_file_stream<<"SNP is not measured for id with number "<<id_num_in_set1<<" for snp with number "<<i<<" in set1. But measured in set2.\n";			
									warnig_file_stream<<"In output set: id_number="<<id_num_in_set1<<", snp_number="<<i<<"\n";
									//std::cout<<"SNP is not measured for id with number "<<id_num_in_set1<<" for snp with number "<<i<<" in set1. But measured in set2.\n";			
									//std::cout<<"In output set: id_number="<<id_num_in_set1<<", snp_number="<<i<<"\n";
	
									if(replace_na == YES)
										{
										Set_merged.set(id_num_in_set1, i, Set2.get(k,snp_in_set2));
										warnig_file_stream<<"Vallue of this field is taken from vallue from same field in set2\n";
										//std::cout<<"Vallue of this field is taken from vallue from same field in set2\n";
										}
									else
										{
										//		std::cout<<"OOOO\n";
										}
									warnig_file_stream<<"\n";
									}
								}
							else
								{
	
								if(int(Set1.get(id_num_in_set1, i)) != 0 &&  int(Set2.get(k, snp_in_set2)) == 0)
									{
									warnig_file_stream<<"SNP is not measured for id with number "<<k<<" for snp with number "<<snp_in_set2<<" in set2. But measured in set1.\n";			
									warnig_file_stream<<"In output set: id_number="<<k - shift_k + num_ids1<<", snp_number="<<i<<"\n";
									//std::cout<<"SNP is not measured for id with number "<<k<<" for snp with number "<<snp_in_set2<<" in set2. But measured in set1.\n";			
									//std::cout<<"In output set: id_number="<<k - shift_k + num_ids1<<", snp_number="<<i<<"\n";
									if(replace_na == YES)
										{
										warnig_file_stream<<"Vallue of this field is taken from vallue from same field in set1\n";
										//std::cout<<"Vallue of this field is taken from vallue from same field in set1\n";
										}
									else
										{
										Set_merged.set(id_num_in_set1, i, 0);
										}
									warnig_file_stream<<"\n";
									}
								else
									{
									Set_merged.set(id_num_in_set1, i, Set2.get(k, snp_in_set2));
									}

							
								}
							}



						shift_k++;
					 	continue;
						}

//					std::cout<<"Set2.get("<<k<<","<<snp_in_set2<<")="<<int(Set2.get(k,snp_in_set2))<<"\n";
//					std::cout<<"k="<<k<<", shift_k="<<shift_k<<", num_ids1="<<num_ids1<<"\n";
//					std::cout<<"Set_merged.set("<<k - shift_k + num_ids1<<", "<<i<<", "<<int(Set2.get(k,snp_in_set2))<<")\n";
					Set_merged.set(k - shift_k + num_ids1, i, Set2.get(k,snp_in_set2));
					}
				}



//			else
//				{	
//				for(int k=num_ids1+1 ; k<=num_ids_total_in_new_set ; k++)
//					{
//					printf(_("12\n"));
//					Set_merged.set(k, i, 0); //Do you realy need it??? Comment it if you are sure that R pass object with all zeros
//					std::cout<<"Set_merged.set("<<k<<", "<<i<<", "<<0<<")"<<"\n";
//					printf(_("13\n"));
//					}
//				}
			//----------------------
			//----------------------





		}
	else
		{
		unsigned shift_k=0;
		
		
//		std::cout<<"13.0: i="<<i<<", shift_i="<<shift_i<<"\n";;
		if(search.is_it_snp_in_set1(i-num_snps1))
	 		{
//			std::cout<<"!!!search.is_it_snp_in_set1("<<i-num_snps1<<")="<<search.is_it_snp_in_set1(i-num_snps1)<<"\n"; 
//			std::cout<<"13.0.1\n";
			shift_i++; 
//			std::cout<<"13.0.2\n";
			continue;
			} // This snp has storaged already for both sets.
		
//		std::cout<<"13.1\n";
			

//		std::cout<<"i="<<i<<"\n";
		for(unsigned k=1 ; k<=num_ids2 ; k++)
			{
//			std::cout<<"k="<<k<<"\n";
//			printf(_("18\n"));
			if(unsigned num_id_in_set1 = search.what_id_is_in_set1(k))
					{
//					std::cout<<"_____________________________\n";
//					std::cout<<"set_mereged[29]="<<int(set_mereged[29])<<"\n";
//					std::cout<<"set_mereged[0]="<<int(set_mereged[0])<<"\n";
//					std::cout<<"search.what_id_is_in_set1("<<k<<")="<<search.what_id_is_in_set1(k)<<"\n";
//					std::cout<<"Set2.get("<<k<<","<< i-num_snps1<<")="<<int(Set2.get(k, i-num_snps1))<<"\n";
//					std::cout<<"Set_merged.set("<<num_id_in_set1<<","<<i-shift_i<<","<<int(Set2.get(k, i-num_snps1))<<")"<<"\n";
//					std::cout<<"Set_merged.get("<<num_id_in_set1<<","<<i-shift_i<<")="<<int(Set_merged.get(num_id_in_set1, i-shift_i))<<"\n";
					Set_merged.set(num_id_in_set1, i-shift_i, Set2.get(k, i-num_snps1));
//					std::cout<<"Set_merged.get("<<num_id_in_set1<<","<<i-shift_i<<")="<<int(Set_merged.get(num_id_in_set1, i-shift_i))<<"\n";
					shift_k++;
//					std::cout<<"_____________________________\n";
					}
			else	
				{
//				std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
//				std::cout<<"Set_merged.set("<<k-shift_k+num_ids1<<","<< i-shift_i<<","<<int(Set2.get(k, i-num_snps1))<<")"<<"\n";	
//				std::cout<<"k="<<k<<", shift_k="<<shift_k<<", num_ids1="<<num_ids1<<", i="<<i<<", shift_i="<<shift_i<<"\n"; 
//					std::cout<<"Set_merged.get("<<k-shift_k+num_ids1<<","<<i-shift_i<<")="<<int(Set_merged.get(k-shift_k+num_ids1, i-shift_i))<<"\n";
				Set_merged.set(k-shift_k+num_ids1, i-shift_i, Set2.get(k, i-num_snps1));
//					std::cout<<"Set_merged.get("<<k-shift_k+num_ids1<<","<<i-shift_i<<")="<<int(Set_merged.get(k-shift_k+num_ids1, i-shift_i))<<"\n";
//				std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
				}
		
			}
			
//		std::cout<<"i="<<i<<", num_snps1="<<num_snps1<<"\n";
		}


//if(i == 8) return;
	
	}//end of for cycle

warnig_file_stream.close();



std::cout<<"End of Program!!!\n";







//		error(_("'ord' must be a positive integer")); //send message to R 
//		printf(_("\nAll right?????????????\n"));
//		Rp33333333333333dgdsgdsgsdg\n"));



	}
//------------------------------------------------------------








}
