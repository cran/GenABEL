#ifndef SMV_GTPS_CONTAINER_H
#define SMV_GTPS_CONTAINER_H














//==================================================================
// Class gtps_container - container for gtps data. Class facilitates access to gtps array.
class gtps_container
	{
	public:
		//gtps_array_raw - raw array is gotten from R, id_numbers - ID amount, snp_numbers - SNP amount
		gtps_container(char * gtps_array_raw, unsigned id_numbers, unsigned snp_numbers);
		~gtps_container(void);

		//Get polymorphism for person who have ID=id_position for SNP=snp_position
		//Here the first id_position equal zerro, the second - one etc...
		char get(unsigned id_position, unsigned snp_position);
		
		
		void set(unsigned id_position, unsigned snp_position, char data);

	private:
		void get_our_byte_number_and_local_person_number(unsigned id_position, unsigned snp_position);

		char *gtps_array; //pointer to array where we stotages our data (passed from R)
		unsigned id_numbers, snp_numbers;
		unsigned nbytes_for_one_snp;

		unsigned our_byte_number, 
						 local_person_number; //can have vallues: 1, 2, 3 or 4. Show what is position in the byte.
//	  const unsigned rearrangement_array[];
		unsigned *rearrangement_array;


	};
//==================================================================

#endif
