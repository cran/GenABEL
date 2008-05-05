//=====================================================================================
//
//     		   Filename:  convert.snp.affymetrix.cpp
//
//		    Description:  Function for converting affymetrix data to GenABEL raw data format
//
//    		    Version:  1.0
//      		  Created:  16-Apr-2008
//    		   Revision:  none
//	last modification: 16-Apr-2008
//
//  		       Author:  Maksim V. Struchalin, Yurii S. Aulchenko
//     		 	  Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
//       		   	Email:  m.struchalin@erasmusmc.nl, i.aoultchenko@erasmusmc.nl
//
//=====================================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include "Chip.h"
//#include "/home/maksim/work/GenABEL_dev/GenABEL/src/gtps_container.h"

#include <R.h>  // to include Rconfig.h

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pkg", String)
// replace pkg as appropriate
#else
#define _(String) (String)
#endif





extern "C" {


std::string replace(std::string val, char what, char replace);



void convert_snp_affymetrix_C(char **dirname_, char **filelist, unsigned *files_amount_, char **map_filename_, char **outfilename_, unsigned *skipaffym, char **alleleID_names, char *alleleID, unsigned *alleleID_amount)
{

char *outfilename = *outfilename_;
char *dirname = *dirname_;
char *map_filename = *map_filename_;


unsigned files_amount=*files_amount_;



std::map<std::string, char> coding;
for(unsigned i=0 ; i<*alleleID_amount ; i++)
	{
	coding[alleleID_names[i]] = alleleID[i];
	}


std::cout<<"reading map...\n";
AffymetrixChipMap Map(map_filename, 2, 0, 2, 4, 5, 3, 9, 10, 6);
std::cout<<"map is read...\n";





std::vector<ChipData *> ids_chip;
for(unsigned i=0 ; i<files_amount ; i++)
	{
	const char * file = (std::string(dirname) + "/" + std::string(filelist[i])).c_str();
	std::cout<<"opening file "<<file<<"...\n";
	ids_chip.push_back(new affymetrix_chip_data(file, 0, 1, *skipaffym));
	}





unsigned id_amount=ids_chip.size(); 


std::ofstream outfile(outfilename);
if(!outfile.is_open()){std::cout<<"Can not open file "<<outfilename<<"\n"; exit(1);}

std::cout<<"Save to file "<<outfilename<<"\n";


outfile << "#GenABEL raw data version 0.1\n";

//save IDs
std::cout<<"saving Id names...\n";
for(unsigned id=0 ; id<files_amount ; id++)
	{
	outfile<<replace(std::string(filelist[id]), ' ', '_')<<" ";
	}
outfile<<"\n";

std::string snpname;

//save snpnames
std::cout<<"saving SNP names...\n";
unsigned snp_amount=ids_chip[0]->get_snp_amount();
for(unsigned snp=0 ; snp<snp_amount ; snp++)
	{
	snpname = ids_chip[0]->get_snp_name(snp);
	if(Map.is_snp_in_map(snpname)){outfile<<Map.recode_snp(snpname.c_str())<<" ";}
	}
outfile<<"\n";

//save chromosome 
std::cout<<"saving chromosome data...\n";
for(unsigned snp=0 ; snp<snp_amount ; snp++)
	{
	snpname = ids_chip[0]->get_snp_name(snp);
	if(Map.is_snp_in_map(snpname)){outfile<<Map.get_chromosome(snpname.c_str())<<" ";}
	}
outfile<<"\n";


//save position (map) 
std::cout<<"saving position data...\n";
for(unsigned snp=0 ; snp<snp_amount ; snp++)
	{
	snpname = ids_chip[0]->get_snp_name(snp);
	if(Map.is_snp_in_map(snpname)){outfile<<Map.get_phisical_position(snpname.c_str())<<" ";}
	}
outfile<<"\n";



//save coding
std::cout<<"saving coding data...\n";
outfile.flags(std::ios_base::hex); //for what is it <-?
for(unsigned snp=0 ; snp<snp_amount ; snp++)
	{
	snpname = ids_chip[0]->get_snp_name(snp);
	if(Map.is_snp_in_map(snpname))
		{
		outfile.width(2);
		outfile.fill('0');
		static std::string allele_A, allele_B;
		allele_A = Map.get_allele_A(snpname.c_str());
  	allele_B = Map.get_allele_B(snpname.c_str());
		outfile<<unsigned(coding[allele_A+allele_B])<<" ";
		}
	}
outfile<<"\n";





//save strand
std::cout<<"saving strand data...\n";
std::map<char, unsigned> strand_recode;
strand_recode['u']=0;
strand_recode['+']=1;
strand_recode['-']=2;

for(unsigned snp=0 ; snp<snp_amount ; snp++)
	{
	snpname = ids_chip[0]->get_snp_name(snp);
	if(Map.is_snp_in_map(snpname))
		{
		outfile.width(2);
		outfile.fill('0');
		static char strand;
		strand = Map.get_strand(snpname.c_str());
		outfile<<strand_recode[strand]<<" ";
		}
	}
outfile<<"\n";



//save polymorphism data
std::cout<<"saving polymorphism data...\n";
unsigned long gtps_byte_amount = (unsigned long)ceil((double)id_amount/4.);
char *gtps_for_one_snp = new char[gtps_byte_amount];



unsigned *rearrangement_array = new unsigned[4];
rearrangement_array[0] = 6;
rearrangement_array[1] = 4;
rearrangement_array[2] = 2;
rearrangement_array[3] = 0;


for(unsigned snp=0 ; snp<snp_amount ; snp++)
	{
	if(!Map.is_snp_in_map(ids_chip[0]->get_snp_name(snp))) {continue;} // skip SNP if it doesn't exsist in our MAP
	for(unsigned i=0 ; i<gtps_byte_amount ; i++) gtps_for_one_snp[i]=0;
	
	static unsigned counter1, counter2;
	counter1=counter2=0;
	
	for(unsigned id=0 ; id<id_amount ; id++)
		{
		gtps_for_one_snp[counter2] = gtps_for_one_snp[counter2] | ids_chip[id]->get_polymorphism(snp)	<< rearrangement_array[counter1];
		counter1++;
		if(counter1==4) {counter1=0; counter2++;}
		}
		

	for(unsigned id=0 ; id<gtps_byte_amount ; id++)
		{
		outfile.width(2);
  	outfile.fill('0');
		outfile<<unsigned(gtps_for_one_snp[id]&0xFF)<<" ";
		}
	outfile<<"\n";
	}

delete gtps_for_one_snp;	
delete rearrangement_array;


if(Map.get_exclude_amount() != 0) std::cout<<Map.get_exclude_amount()<<" SNPs excluded because of absent enough information in map file\n";

std::cout<<"Finshed... Data saved into file "<<outfilename<<"\n";
outfile.close();
}











//_________________________________________________________________
std::string replace(std::string val, char what, char replace='_')
{
unsigned length = val.length();

for(unsigned i=0 ; i<length ; i++)
	{
	if(val[i]==what) val[i]=replace;
	}
return val;
}
//_________________________________________________________________

}//end of externl C
