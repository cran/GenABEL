/**
 *
 * 2007.07.26 by Yurii Aulchenko, EMCR
 *
 * last modified 2007.07.26 by YA
 *
 **/

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <stdio.h>

using namespace std;

//STL
#include <vector>
#include <vector>
#include <map>
#include <iterator>

#include <R.h>

#define MAXSNPs 2000000

/**
 *
 * strand == 0 -> unknown ('u')
 * strand == 1 -> plus ('+')
 * strand == 2 -> minus ('-')
 * strand == 3 -> four columns (name,chr,pos,strand) expected at the beginning
 *
 **/

extern "C" {
void convert_snp_illumina (char** filename, char** outfilename, int* Strandid, int* Bcast, char **allele_codes, int* Ncodes) {

	short unsigned int ncodes = *Ncodes;
	short unsigned int strandid = *Strandid;
	long unsigned int bcast = *Bcast;

	bool loud = TRUE;
	//if (bcast <= 0) loud = FALSE;

	long unsigned int nsnps=0;
	long unsigned int nids=0;
	unsigned long int nbytes=0;
	unsigned long int byte;

	string data;
	string tempstr;

	vector<string> iid;
	vector<string> chrom; string tmp_chrom;
	vector<string> snpnam; string tmp_snpnam;
	vector<unsigned long> map; unsigned long tmp_map;
	vector<unsigned char*> gtype; unsigned char* tmp_gtype;
	vector<string> coding; string tmp_coding,tmp_coding1;
	vector<unsigned short int> intcoding;
	vector<unsigned short int> strand; string tmp_strand;
	vector<string> codeset(ncodes);
	char tmp_chcoding [10];

	for (int i=0;i<ncodes;i++) codeset[i].assign(allele_codes[i]);

	ifstream illfile (filename[0]);
	if (!illfile.is_open()) {
		error ("could not open file '%s'!",filename[0]);
	}

	if (loud) {
		Rprintf("Reading genotypes from file '%s' ...\n",filename[0]);
	}

	nsnps = 0;
	unsigned long int lasti = 1;

	//
	// Processing header line
	//

	if (getline(illfile,data)) {
		istringstream datas (data);
		datas >> tempstr >> tempstr >> tempstr;
		if (strandid == 3) datas >> tempstr;
		while (datas >> tempstr) {
			iid.push_back(tempstr);
			nids++;
		}
	} else {
		error("Can not read the first line from file '%s'!\n",filename[0]);
	}

	char* chgt = new char [nids*2];
	nbytes = (unsigned long int)ceil((double)nids/4.);

	//
	//Processing rest of the file
	//
		
	char gdata;

	while (getline(illfile,data)) {
		nsnps++;

		istringstream datas (data);


		if (!(datas >> tmp_snpnam >> tmp_chrom >> tmp_map))
			error("First three fields are missing in line %i!\n",(nsnps+1));
		else {
			chrom.push_back(tmp_chrom);
			snpnam.push_back(tmp_snpnam);
			map.push_back(tmp_map);
			if (strandid == 3) {
				if (!(datas >> tmp_strand)) error("Strand field is missing in line %i, SNP '%s'!\n",(nsnps+1),tmp_snpnam.c_str());
				if (tmp_strand == "u") strand.push_back(0);
				else if (tmp_strand == "+") strand.push_back(1);
				else if (tmp_strand == "-") strand.push_back(2);
				else error("Bad strand coding ('%s'), only '+', '-' or 'u' is accepted!\n",tmp_strand.c_str());
			} else strand.push_back(strandid);

			unsigned long int idx = 0;

			idx = 0;
			for (unsigned long int i = 0;i<nids;i++) {
				if (!(datas >> tempstr)) error("Too few fields for SNP '%s', line %i (no. IDs = %li, no. fields = %li)!\n",snpnam[nsnps-1].c_str(),(nsnps+1),nids,(i+1));
				char sd[10];
				sprintf(sd,"%s",tempstr.c_str());
				chgt[idx++] = sd[0];
				chgt[idx++] = sd[1];
				//std::cout<<"sd[0]="<<sd[0]<<", sd[1]="<<sd[1]<<"\n";
			}

			//
			//Processing SNP string
			//
			unsigned short int* gnum = new unsigned short int [nids*2];
			unsigned short int ind;
			unsigned short int offset[4] = {6,4,2,0};
			char allele1 = 0;
			char allele2 = 0;
			unsigned long int ca1 = 0;
			unsigned long int ca2 = 0;


			for (unsigned long int idx = 0; idx < 2*nids; idx++) 
				{
				gdata = chgt[idx];
				if (gdata == allele1) 
					{
					gnum[idx] = 1;
					ca1++;
					}
			 	else if (gdata == allele2) 
					{
					gnum[idx] = 3;
					ca2++;
					}
			 	else if (gdata == '0' || gdata == '-') 
					{
					gnum[idx] = 0;
					}
			 	else 
					{
					if (allele1 == 0)
				 		{
						allele1 = gdata;
						gnum[idx] = 1;
						ca1++;
						}
				 	else if (allele2 == 0) 
						{
						allele2 = gdata;
						gnum[idx] = 3;
						ca2++;
						}
				 	else 
						{
						error ("illegal genotype (three alleles) for SNP '%s' (line %li)!", snpnam[nsnps-1].c_str(),(nsnps+1));
						}
					}
				//std::cout<<"gnum["<<idx<<"]="<<gnum[idx]<<"\n";
				}

			if (!allele1 && !allele2) tmp_coding="12"; // all genotypes missing
			else if (!allele1 && allele2) sprintf(tmp_chcoding,"%c%c",allele2,allele2); // only one allele present
			else if (allele1 && !allele2) sprintf(tmp_chcoding,"%c%c",allele1,allele1); // only one allele present
			else if (allele1 && allele2) {
				if (ca1 > ca2) sprintf(tmp_chcoding,"%c%c",allele1,allele2);
				else sprintf(tmp_chcoding,"%c%c",allele2,allele1);
			}
			//Rprintf("%s\n",tmp_chcoding);
			tmp_coding.assign(tmp_chcoding);
			int ccd = -1;
			for (int i = 0; i < ncodes; i++) {
				if (codeset[i].compare(tmp_coding)==0) {
					ccd = i + 1;
					intcoding.push_back(ccd);
				}
			}
			if (ccd<0) error ("coding '%s' for SNP not recognised !\n",tmp_coding.c_str());
			try {
				tmp_gtype = new unsigned char [nbytes];
			}
			catch (bad_alloc) {
				error ("ran out of memory reading file '%s' line %li !");
			}
			if (tmp_gtype == NULL) {
				error ("ran out of memory reading file '%s' line %li !");
			}

			idx = 0;
			for (byte = 0; byte < nbytes; ++byte) {

				tmp_gtype[byte] = 0;
				for (ind = 0; ind < 4; ++ind) {

					switch (gnum[idx]+gnum[idx+1]) {
					case 2:
						if (ca1 > ca2)
							tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)1 << offset[ind]);
						else
							tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)3 << offset[ind]);
						break;
					case 4:
						tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)2 << offset[ind]);
						break;
					case 6:
						if (ca1 > ca2)
							tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)3 << offset[ind]);
						else
							tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)1 << offset[ind]);
						break;
					case 0:
						tmp_gtype[byte] = tmp_gtype[byte] | (0 << offset[ind]); // this does nothing
						break;
					default:
						error ("illegal genotype (half missing) SNP '%s' file '%s' line %li !",
								snpnam[nsnps-1].c_str(),filename[0],(nsnps-1));
					}
					idx += 2;
					if (idx >= 2*nids) break;
				}
			
			//std::cout<<"tmp_gtype["<<byte<<"]="<<int(tmp_gtype[byte])<<"\n";
			}

			gtype.push_back(tmp_gtype);

			if (loud && ((nsnps)*nids) > ((bcast)*lasti)) {
				Rprintf("  ... analysed %li genotypes ...\n",((nsnps)*nids));
				lasti += 1;
			}

		delete [] gnum;
//		delete [] tmp_gtype;   We should not delete it here because it stores our data. Delete it at the very end.
		}
	
	
	
	}

	const ios_base::fmtflags hex = ios_base::hex;

	ofstream outfile (outfilename[0]);
	if (!outfile.is_open()) {
		error ("could not open file '%s' !",outfilename[0]);
	}

	if (loud) {
		Rprintf("Writing to file '%s' ...\n",outfilename[0]);
	}

	outfile << "#GenABEL raw data version 0.1";
	outfile << endl;

	copy(iid.begin(), iid.end(), ostream_iterator<string>(outfile, " "));
	outfile << endl;

	copy(snpnam.begin(), snpnam.end(), ostream_iterator<string>(outfile, " "));
	outfile << endl;

	copy(chrom.begin(), chrom.end(), ostream_iterator<string>(outfile, " "));
	outfile << endl;

	copy(map.begin(), map.end(), ostream_iterator<unsigned long>(outfile, " "));
	outfile << endl;

	outfile.flags(hex);

	for (unsigned long int i=0;i<nsnps;i++) {
		outfile.width(2);
		outfile.fill('0');
		outfile << (unsigned int)intcoding[i] << " ";
	}
	outfile << endl;

	for (unsigned long int i=0;i<nsnps;i++) {
		outfile.width(2);
		outfile.fill('0');
		outfile << (unsigned int)strand[i] << " ";
	}
	outfile << endl;

	for (unsigned long int i=0;i<nsnps;i++) {

		for (byte = 0; byte < nbytes; ++byte) {
			outfile.width(2);
			outfile.fill('0');
			outfile << (unsigned int)gtype[i][byte];
			//std::cout<<"(unsigned int)gtype["<<i<<"]["<<byte<<"]="<<(unsigned int)gtype[i][byte]<<" ";
			outfile << " ";
		}
		outfile << endl;
		//std::cout<<"\n";

	}

	if (loud) {
		Rprintf("... done.\n");
	}

	for (unsigned long int i=0;i<nsnps;i++) 
		{
		delete [] gtype[i];
		}


    delete [] chgt;
}
}
