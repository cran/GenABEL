#include "gtps_container.h"
#include <cstdlib>  
#include <math.h>













gtps_container::gtps_container(char * gtps_array_, unsigned id_numbers_, unsigned snp_numbers_)
{
rearrangement_array = new unsigned[4];
rearrangement_array[0] = 6;
rearrangement_array[1] = 4;
rearrangement_array[2] = 2;  
rearrangement_array[3] = 0;



gtps_array = gtps_array_;	
id_numbers = id_numbers_;
snp_numbers = snp_numbers_;
nbytes_for_one_snp = unsigned(ceil(double(id_numbers)/4.) + 0.5);

//std::cout<<"gtps_container::gtps_container: nbytes_for_one_snp="<<nbytes_for_one_snp<<"\n";

}
//------------------------------------------------------------------


gtps_container::~gtps_container(void)
{
delete[] rearrangement_array;
}




//------------------------------------------------------------------
void gtps_container::get_our_byte_number_and_local_person_number(unsigned id_position, unsigned snp_position)
{
//std::cout<<"\n---------------------------------\n";
//std::cout<<"id_position="<<id_position<<", snp_position="<<snp_position <<"\n";
//std::cout<<"byte_amount_for_one_snp="<<nbytes_for_one_snp<<", id_numbers="<<id_numbers<<"\n";	
//std::cout<<" int(ceil(id_position/4.)+0.5)="<<int(ceil(id_position/4.)+0.5)<<"\n";
//std::cout<<"(snp_position-1)*byte_amount_for_one_snp="<<(snp_position-1)*nbytes_for_one_snp<<"\n";
our_byte_number = int(ceil(id_position/4.)+0.5) + (snp_position-1)*nbytes_for_one_snp;   //What byte is our id in? 

//std::cout<<"our_byte_number="<<our_byte_number<<"\n";
//std::cout<<"---------------------------------\n";

local_person_number = id_position - ((our_byte_number-(snp_position-1)*nbytes_for_one_snp-1)*4); //What position of our id is in byte?
}
//------------------------------------------------------------------









//------------------------------------------------------------------
char gtps_container::get(unsigned id_position, unsigned snp_position)
{
//static const unsigned rearrangement_array[] = {6, 4, 2,0}; 

//std::cout<<"gtps_container::get: START:\n";
//std::cout<<"gtps_container::get: id_position="<<id_position<<", snp_position="<<snp_position<<"\n";
get_our_byte_number_and_local_person_number(id_position, snp_position);
//std::cout<<"gtps_container::get: our_byte_number="<<our_byte_number<<", local_person_number="<<local_person_number<<"\n";




char our_byte_vallue = gtps_array[our_byte_number-1];


//std::cout<<"gtps_container::get: our_byte_vallue="<<int(our_byte_vallue)<<", rearrangement_array[local_person_number-1])="<<rearrangement_array[local_person_number-1]<<", local_person_number="<<local_person_number<<"\n";





//for(int i=0 ; i<4 ; i++)
//std::cout<<"rearrangement_array["<<i<<"]="<<rearrangement_array[i]<<"\n";




return char((gtps_array[our_byte_number-1] >> rearrangement_array[local_person_number-1]) & 3); 
}
//------------------------------------------------------------------











//------------------------------------------------------------------
void gtps_container::set(unsigned id_position, unsigned snp_position, char data)
{
//for(int i=0 ; i<4 ; i++)
//std::cout<<"rearrangement_array["<<i<<"]="<<rearrangement_array[i]<<"\n";

		
static const char clear_info_for_person[]={63, 207, 243, 252};
		
//std::cout<<"set: id_position="<<id_position<<", snp_position="<<snp_position<<", data="<<int(data)<<"\n";
get_our_byte_number_and_local_person_number(id_position, snp_position);


//std::cout<<"set: (gtps_array[our_byte_number]&clear_info_for_person[local_person_number-1])="<<int(gtps_array[our_byte_number]&clear_info_for_person[local_person_number-1])<<"\n";

//char tmp = (gtps_array[our_byte_number]&clear_info_for_person[local_person_number-1]) | data;
//std::cout<<"set: tmp="<<int(tmp)<<"\n";


//std::cout<<"gtps_array["<<our_byte_number-1<<"] = (gtps_array["<<our_byte_number-1<<"]&"<<int(clear_info_for_person[local_person_number-1])<<" |"<< int((data&3)<<rearrangement_array[local_person_number-1])<<"\n";

gtps_array[our_byte_number-1] = (gtps_array[our_byte_number-1]&clear_info_for_person[local_person_number-1]) | (data&3)<<rearrangement_array[local_person_number-1];
}
//------------------------------------------------------------------




/*
gtps_container::checking(unsigned id_position, unsigned snp_position)
{
if (id_position > id_numbers || id_position<0)
	{
	std::cerr<<"Bad ID number. ID_number="<<id_position<<". Maximum ID_number can be "<<id_position<<"\n";
	}

if(snp_position > snp_position || snp_position <0)
	{
	std::cerr<<"Bad ID number. ID_number="<<id_position<<". Maximum ID_number can be "<<id_position<<"\n";
	}
}
*/



