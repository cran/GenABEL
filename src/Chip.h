#ifndef SMV_CHIP
#define SMV_CHIP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>


class ChipData
	{
	public:
//		virtual ~ChipData()=0;
		virtual const unsigned get_snp_amount(void)=0;
		virtual const int get_polymorphism(unsigned snp_num)=0; 
		virtual const char * get_snp_name(unsigned snp_num)=0; 
	
	
	};



class affymetrix_chip_data : public ChipData
	{
	public:
		affymetrix_chip_data(std::string filename, unsigned snp_positionm, unsigned polymorphism_position,  unsigned skip_first_line);
		~affymetrix_chip_data();
	
		const unsigned get_snp_amount(void);
		const int get_polymorphism(unsigned snp_num); //0->NA, 1->AA, 2->AB, 3->BB 
		const char * get_snp_name(unsigned snp_num); //rsp12345678


	private:
		std::string filename;
		unsigned snp_amount;
//		std::vector<char> polymorphism; //0,1,2,3
		char * polymorphism; //0,1,2,3
//		std::vector<std::string> snp_name; 			
		char ** snp_name; 			
				

	};





//__________________________________________________________
struct map_values
	{
	std::string snp_name; //affymetrix SNP name
	std::string recoded_snp_name; //rs12345678
	std::string phisical_position; //12345678
	char strand; //+ or - or u or another
	std::string chromosome;
	std::string	allele_A; //A, T, G, C
	std::string allele_B; //A, T, G, C
	};


class ChipMap
	{
	public:
//		virtual ~ChipMap()=0;
	
		virtual std::string recode_snp(const char*); //recode "SNP_A-1780619" -> "rs17106009"
		virtual std::string get_phisical_position(const char*); 
		virtual char get_strand(const char*);
		virtual std::string get_chromosome(const char*);
		virtual std::string get_allele_A(const char*);
		virtual std::string get_allele_B(const char*);
		virtual bool is_snp_in_map(std::string);
	protected:
		std::map<std::string, map_values> Map;
		

	};






class AffymetrixChipMap : public ChipMap
	{
	public:		
		AffymetrixChipMap(const char* filename, unsigned skip_first_lines, unsigned snp_name_position, unsigned recoded_snp_name_position, unsigned phisical_position_position, unsigned strand_position, unsigned chromosome_position, unsigned allele_A_position, unsigned allele_B_position, unsigned reg1_position, char delim=',');
		~AffymetrixChipMap();
		unsigned get_exclude_amount(void);
	private:
		unsigned exclude_amount;
	};

std::string cut_quotes(std::string str);




#endif //SMV_CHIP 
