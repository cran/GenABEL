/**
 *
 * Code contributed by Toby Johnson
 *
 * modified by Yurii Aulchenko:
 * list -> vector
 * output hex format -- ensure that two chars are output for every number (width=2)
 * conversion to GenABEL raw format version 0.1 
 *
 * last modified 2007.12.18
 *
 **/
#include <cstdlib>  
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std; 

//STL
#include <vector>
#include <vector>
#include <map>
#include <iterator>

#include <R.h>




extern "C" {
void convert_snp_tped (char** tpedfilename, char** tfamfilename, char** outfilename, int* Strandid, int* bcast, char **allele_codes, int* Ncodes) {

	short unsigned int ncodes = *Ncodes;
	short unsigned int strandid = *Strandid;

	int verbose = *bcast ? 1 : 0;

	long int linecount = 0;
	string data;
	string token;

	vector<string> iid; string tmp_iid;
	vector<string> chrom; string tmp_chrom;
	vector<string> snpnm; string tmp_snpnm;
	vector<double> genmap; double tmp_genmap;
	vector<unsigned long> phymap; unsigned long tmp_phymap;
	vector<unsigned char*> gtype; unsigned char* tmp_gtype;
	vector<string> coding; string tmp_coding,tmp_coding1;
	vector<unsigned short int> intcoding;
	vector<unsigned short int> strand; string tmp_strand;
	vector<string> codeset(ncodes);
	char tmp_chcoding [10];

	for (int i=0;i<ncodes;i++) codeset[i].assign(allele_codes[i]);

	///////////////////////
	// read the tfamfile //
	///////////////////////

	ifstream tfamfile (tfamfilename[0]);
	if (!tfamfile.is_open()) {
		error ("could not open file '%s' !",tfamfilename[0]);
	}

	if (verbose) {
		Rprintf("Reading individual ids from file '%s' ...\n",tfamfilename[0]);
	}

	while (getline(tfamfile,data)) {
		istringstream datas (data);
		if ( datas >> token >> tmp_iid) {
			iid.push_back(tmp_iid);
		}
	}

	tfamfile.close();

	int nids = iid.size();
	int nbytes = (int)ceil((double)nids/4.);

	if (verbose) {
		Rprintf("... done.  Read %i individual ids from file '%s'\n",nids,tfamfilename[0]);
	}

	///////////////////////
	// read the tpedfile //
	///////////////////////

	int idx;
	char gdata;
	int* gnum = new int [nids*2];

	int byte;
	int ind;
	int offset[4] = {6,4,2,0};


	ifstream tpedfile (tpedfilename[0]);
	if( !tpedfile.is_open()) {
		error ("could not open file '%s' !",tpedfilename[0]);
	}

	if (verbose) {
		Rprintf("Reading genotypes from file '%s' ...\n",tpedfilename[0]);
	}

	linecount = 0;
	while (getline(tpedfile,data)) {
		++linecount;

		istringstream datas (data);

		if (datas >> tmp_chrom >> tmp_snpnm >> tmp_genmap >> tmp_phymap) {
			chrom.push_back(tmp_chrom);
			snpnm.push_back(tmp_snpnm);
			genmap.push_back(tmp_genmap);
			phymap.push_back(tmp_phymap);
			strand.push_back(strandid);

			char allele1 = 0;
			char allele2 = 0;
			unsigned long int ca1 = 0;
			unsigned long int ca2 = 0;

			for (idx = 0; idx < 2*nids; ++idx) {
				if (datas >> gdata) {

					if (gdata == allele1) {
						gnum[idx] = 1;
						ca1++;
					} else if (gdata == allele2) {
						gnum[idx] = 3;
						ca2++;
					} else if (gdata == '0') {
						gnum[idx] = 0;
					} else {
						if (allele1 == 0) {
							allele1 = gdata;
							gnum[idx] = 1;
						} else if (allele2 == 0) {
							allele2 = gdata;
							gnum[idx] = 3;
						} else {
							error ("illegal genotype (three alleles) snp '%s' file '%s' line %li !",
									tmp_snpnm.c_str(),tpedfilename[0],linecount);
						}
					}

				} else {
					error ("not enough genotypes for snp '%s' file '%s' line %li !",
							tmp_snpnm.c_str(),tpedfilename[0],linecount);
				}
			}


			//	if (ca1 > ca2) sprintf(tmp_chcoding,"%c%c",allele1,allele2);
			//	else sprintf(tmp_chcoding,"%c%c",allele2,allele1);
			if (!allele1 && !allele2) tmp_coding="12"; // all genotypes missing
			else if (!allele1 && allele2) sprintf(tmp_chcoding,"%c%c",allele2,allele2); // only one allele present
			else if (allele1 && !allele2) sprintf(tmp_chcoding,"%c%c",allele1,allele1); // only one allele present
			else if (allele1 && allele2) {
				if (ca1 > ca2) sprintf(tmp_chcoding,"%c%c",allele1,allele2);
				else sprintf(tmp_chcoding,"%c%c",allele2,allele1);
			}
			//Rprintf("%s\n",tmp_chcoding);
			tmp_coding.assign(tmp_chcoding);
			//	if (!allele1 || !allele2) tmp_coding="12";
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
								snpnm[linecount-1].c_str(),tpedfilename[0],(linecount-1));
					}
					idx += 2;
					if (idx >= 2*nids) break;
				}
			}

			gtype.push_back(tmp_gtype);

			if (verbose && linecount % *bcast == 0) {
				Rprintf("  ... read %li lines ...\n",linecount);
			}

		} else {
			// not even four fields on the line; raise error?
		}

	}

	tpedfile.close();

	if (verbose) {
		Rprintf("...done.  Read %i SNPs from file '%s'\n",chrom.size(),tpedfilename[0]);
	}

	const ios_base::fmtflags hex = ios_base::hex;

	ofstream outfile (outfilename[0]);
	if( !outfile.is_open()) {
		error ("could not open file '%s' !",outfilename[0]);
	}

	if (verbose) {
		Rprintf("Writing to file '%s' ...\n",outfilename[0]);
	}

	outfile << "#GenABEL raw data version 0.1";
	outfile << endl;

	copy(iid.begin(), iid.end(), ostream_iterator<string>(outfile, " "));
	outfile << endl;

	copy(snpnm.begin(), snpnm.end(), ostream_iterator<string>(outfile, " "));
	outfile << endl;

	copy(chrom.begin(), chrom.end(), ostream_iterator<string>(outfile, " "));
	outfile << endl;

	copy(phymap.begin(), phymap.end(), ostream_iterator<unsigned long>(outfile, " "));
	outfile << endl;

	outfile.flags(hex);

	for (unsigned long int i=0;i<chrom.size();i++) {
		outfile.width(2);
		outfile.fill('0');
		outfile << (unsigned int)intcoding[i] << " ";
	}
	outfile << endl;

	for (unsigned long int i=0;i<chrom.size();i++) {
		outfile.width(2);
		outfile.fill('0');
		outfile << (unsigned int)strand[i] << " ";
	}
	outfile << endl;

	for (unsigned long int i = 0;i<gtype.size();i++)
	{
		tmp_gtype = gtype[i];

		for (byte = 0; byte < nbytes; ++byte) {
			outfile.width(2);
			outfile.fill('0');
			outfile << (unsigned int)tmp_gtype[byte];
			outfile << " ";
		}
		outfile << endl;


		//      gtype.pop_front();
	} //while (!gtype.empty());
		
	delete [] tmp_gtype;
	delete [] gnum;

	if (verbose) {
		Rprintf("... done.\n");
	}

}
}



