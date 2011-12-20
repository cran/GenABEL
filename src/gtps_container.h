//#=====================================================================================
//#
//#       Filename:  gtps_container.h
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
//#          Email:  m.struchalin@@erasmusmc.nl, i.aoultchenko@erasmusmc.nl
//#
//#=====================================================================================

#ifndef SMV_GTPS_CONTAINER_H
#define SMV_GTPS_CONTAINER_H

#include <R.h>

#include<iostream>
#include<vector>
#include<memory>








//==================================================================
// Class gtps_container - container for gtps data. Class facilitates access to gtps array.
class gtps_container
	{
	public:
		gtps_container(char * gtps_array_raw, char * strand_array, char *  coding_array, unsigned id_numbers, unsigned snp_numbers);
		gtps_container(char * gtps_array_raw, unsigned id_numbers, unsigned snp_numbers);
		~gtps_container(void);

		//Get polymorphism for person who have ID=id_position for SNP=snp_position
		//Here the first id_position equal zerro, the second - one etc...
		//id_position=1 - THE FIRST ARGUMENT (id_position!=O NEVER!!!). snp_position - same
		char get(unsigned id_position, unsigned snp_position);
		char* get_gtps_array_for_snp(unsigned snp_position);
		char get_strand(unsigned snp_position);
		char get_coding(unsigned snp_position);
		
		
		//id_position=1 - THE FIRST ARGUMENT (id_position!=O NEVER!!!). snp_position - same
		void set(unsigned id_position, unsigned snp_position, char data);

	private:
		void get_our_byte_number_and_local_person_number(unsigned id_position, unsigned snp_position);
		bool do_we_have_strand_and_codding_arrays;
		char *gtps_array; //pointer to array where we stotages our data (passed from R)
		char * strand_array;
		char * coding_array;
		unsigned id_numbers, snp_numbers;
		unsigned nbytes_for_one_snp;

		unsigned our_byte_number, 
						 local_person_number; //can have vallues: 1, 2, 3 or 4. Show what is position in the byte.
//	  const unsigned rearrangement_array[];
		unsigned *rearrangement_array;


	};
//==================================================================

#endif
