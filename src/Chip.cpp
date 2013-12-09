//=====================================================================================
//
//     		   Filename:  Chip.cpp
//
//		    Description:  Function for converting affymetrix data to GenABEL raw data format
//
//    		    Version:  1.0
//      		  Created:  16-Apr-2008
//    		   Revision:  none
//	last modification:  Wed May 14 16:39:33 CEST 2008
//
//  		       Author:  Maksim V. Struchalin, Yurii S. Aulchenko
//     		 	  Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
//       		   	Email:  m.struchalin@erasmusmc.nl, i.aoultchenko@erasmusmc.nl
//
//=====================================================================================

#include "Chip.h"
#include <R.h>


//__________________________________________________________
affymetrix_chip_data::affymetrix_chip_data(std::string filename_, unsigned snp_position, unsigned polymorphism_position, unsigned skip_first_lines)
{
unsigned BUFSIZE=10000000;
char * buf = new char[BUFSIZE];

std::vector<char> polymorphism_vec;
std::vector<std::string> snp_name_vec;


filename = filename_;

std::ifstream file(filename.c_str());
if(!file.is_open()){error("Can not open file \"%s\"\n", filename.c_str());}

for(unsigned i=0 ; i<skip_first_lines ; i++) file.getline(buf,BUFSIZE-1); //skip first line

std::string val;


while(!file.eof())
	{
	file.getline(buf,BUFSIZE-1);
	std::stringstream line_stream(buf);

	for(unsigned col=0; line_stream>>val ; col++ )
		{
		if(col==snp_position) {snp_name_vec.push_back(val);}
		if(col==polymorphism_position) 
			{
			if(std::string(val) == std::string("AA") || val == std::string("1")) polymorphism_vec.push_back(1);
			else if(val == std::string("AB") || val == std::string("2")) polymorphism_vec.push_back(2);
			else if(val == std::string("BB") || val == std::string("3") ) polymorphism_vec.push_back(3);
			else polymorphism_vec.push_back(0);
			}
		if(col>=snp_position && col>=polymorphism_position) break;
		}
	}

static unsigned polym_size;
polym_size = polymorphism_vec.size();
polymorphism = new char[polym_size];
for(unsigned i=0 ; i<polym_size ; i++) polymorphism[i]=polymorphism_vec[i];

snp_amount = snp_name_vec.size();
snp_name = new char*[snp_amount];
for(unsigned i=0 ; i<snp_amount ; i++) 
	{
	snp_name[i] = new char[snp_name_vec[i].length()+1];	
	snp_name_vec[i].copy(snp_name[i], std::string::npos);
	snp_name[i][snp_name_vec[i].length()]=0;
	}
		



file.close();

delete[] buf;

}
//__________________________________________________________


//__________________________________________________________
affymetrix_chip_data::~affymetrix_chip_data(void)
{
delete polymorphism;
for(unsigned i=0 ; i<snp_amount ; i++) 
	{
	delete snp_name[i];
	}
delete [] snp_name;
}
//__________________________________________________________


//__________________________________________________________
const unsigned affymetrix_chip_data::get_snp_amount(void)
{
return snp_amount;
}
//__________________________________________________________




//__________________________________________________________
const int affymetrix_chip_data::get_polymorphism(unsigned snp_num)
{
if(snp_num >= snp_amount) {error("file %s: input SNP amount %i is too big. Maximum is %i\n", filename.c_str(), snp_num, snp_amount-1);} 
return int(polymorphism[snp_num]);
}
//__________________________________________________________


//__________________________________________________________
const char * affymetrix_chip_data::get_snp_name(unsigned snp_num_)
{
if(snp_num_ >= snp_amount) {error("file %s: input SNP amount %i is too big. Maximum is %i\n", filename.c_str(), snp_num_, snp_amount-1);} 
return snp_name[snp_num_];
}
//__________________________________________________________













//__________________________________________________________
std::string ChipMap::recode_snp(const char* snp_name)
{
return Map.find(snp_name)->second.recoded_snp_name;
}
//__________________________________________________________


//__________________________________________________________
std::string ChipMap::get_phisical_position(const char* snp_name)
{
return Map.find(snp_name)->second.phisical_position;
}
//__________________________________________________________



//__________________________________________________________
char ChipMap::get_strand(const char* snp_name)
{
return Map.find(snp_name)->second.strand;
}
//__________________________________________________________



//__________________________________________________________
std::string ChipMap::get_chromosome(const char* snp_name)
{
return Map.find(snp_name)->second.chromosome;
}
//__________________________________________________________



//__________________________________________________________
std::string ChipMap::get_allele_A(const char* snp_name)
{
return Map.find(snp_name)->second.allele_A;
}
//__________________________________________________________



//__________________________________________________________
std::string ChipMap::get_allele_B(const char* snp_name)
{
return Map.find(snp_name)->second.allele_B;
}
//__________________________________________________________





//__________________________________________________________
bool ChipMap::is_snp_in_map(std::string snp_name)
{
static std::map<std::string, map_values>::iterator iter_map;
iter_map=Map.find(snp_name.c_str());
if(iter_map==Map.end()) return false;
else return true;
}
//__________________________________________________________






//__________________________________________________________

AffymetrixChipMap::AffymetrixChipMap(const char* filename, unsigned skip_first_lines, unsigned snp_name_position, unsigned recoded_snp_name_position, unsigned phisical_position_position, unsigned strand_position, unsigned chromosome_position, unsigned allele_A_position, unsigned allele_B_position, unsigned reg1_position, char delim)
{
unsigned BUFSIZE=10000000;
char * buf = new char[BUFSIZE];

std::ifstream file(filename);
if(!file.is_open()){error("Can not open file %s\n", filename);}

for(unsigned i=0 ; i<skip_first_lines ; i++) file.getline(buf,BUFSIZE-1); //skip first line


exclude_amount=0;

map_values map_val;

bool XY_yes, exclude;



while(1)
	{
  file.getline(buf,BUFSIZE-1);
  std::stringstream line_stream(buf);
	
	if(file.eof()) break;
	
	XY_yes = false;
	exclude = false;

	for(unsigned col=0; !line_stream.eof() ; col++ )
		{
		line_stream.getline(buf, BUFSIZE-1, delim);
		
		if(col==recoded_snp_name_position)
	 		{
			map_val.recoded_snp_name=cut_quotes(std::string(buf));
		 	if(map_val.recoded_snp_name == "---") {exclude=true; break;}
			}
		else if(col==snp_name_position)
	 		{
			map_val.snp_name=cut_quotes(std::string(buf));
			}
		else if(col==phisical_position_position)
	 		{
			map_val.phisical_position=cut_quotes(std::string(buf));
			if(map_val.phisical_position == "---" || map_val.phisical_position == "NA") {exclude=true; break;}
			}
		else if(col==strand_position) {map_val.strand=cut_quotes(buf)[0];}
		else if(col==chromosome_position) {map_val.chromosome=cut_quotes(std::string(buf));}
		else if(col==allele_A_position)
	 		{
			map_val.allele_A=cut_quotes(std::string(buf));
			if(map_val.allele_A == "-" || map_val.allele_A == "NA") {exclude=true; break;}
			}
		else if(col==allele_B_position) 
			{
			map_val.allele_B=cut_quotes(std::string(buf));
			if(map_val.allele_A == "-" || map_val.allele_A == "NA") {exclude=true; break;}
			}
		
		if(col == reg1_position) if(cut_quotes(std::string(buf))!=std::string("0")) XY_yes = true;
		}
	if(XY_yes)  map_val.chromosome="XY";
	if(!exclude) {Map[map_val.snp_name]=map_val;}
  else {exclude_amount++;}
	}



file.close();
delete[] buf;
}
//__________________________________________________________


//__________________________________________________________
AffymetrixChipMap::~AffymetrixChipMap()
{
}
//__________________________________________________________



//__________________________________________________________
unsigned AffymetrixChipMap::get_exclude_amount(void)
{
return exclude_amount;
}
//__________________________________________________________




//__________________________________________________________
std::string cut_quotes(std::string str)
{
std::string b;

for(unsigned i=0 ; i<str.length()-1 ; i++)
	{
	if(str[i] == '\"') continue;
	b.push_back(str[i]);
	}

return b;
}
//__________________________________________________________
