/**
 *
 * 2007.07.20 by Yurii Aulchenko, EMCR
 * based on code conributed by Toby Johnson (convert_snp_tped.cpp)
 *
 * modified 2007.07.23 by YA
 * last modified 2008.01.10 by YA
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
#include <map>
#include <iterator>

#include <R.h>

#define MAXIDS 200000

string replace_mrl(string in);
string replace_mach(string in);

/**
 *
 * strand == 0 -> unknown ('u')
 * strand == 1 -> plus ('+')
 * strand == 2 -> minus ('-')
 * strand == 3 -> four columns (name,chr,pos,strand) expected at the beginning
 *
 **/

extern "C" {
void convert_snp_merlin (char** pedfilename, char** mapfilename, char** outfilename,
		int* Strandid, int* bcast, char **allele_codes, int* Ncodes,
		int *Fmt, int *Tra, int *MapHasHeaderLine) {

	int verbose = *bcast ? 1 : 0;

	int ncodes = *Ncodes;
	int strandid = *Strandid;
	int format = *Fmt;
	int traits = *Tra;
	int mapHasHeaderLine = *MapHasHeaderLine;

	long int linecount=0;
	string data;
	string token;
				
	char gdata;

	vector<string> iid; string tmp_iid;
	vector<string> chrom; string tmp_chrom;
	vector<string> snpnm; string tmp_snpnm;
	//    vector<double> genmap; double tmp_genmap;
	vector<unsigned long> phymap; unsigned long tmp_phymap;
	vector<unsigned char*> gtype; unsigned char* tmp_gtype;
	vector<string> coding; string tmp_coding,tmp_coding1;
	vector<int> intcoding;
	vector<int> strand; string tmp_strand;
	vector<string> codeset(ncodes);
	char tmp_chcoding [10];
	//    char tmp_chcoding1 [10];

	for (int i=0;i<ncodes;i++) codeset[i].assign(allele_codes[i]);

	///////////////////////
	// read the mapfile //
	///////////////////////

	ifstream mapfile (mapfilename[0]);
	if( !mapfile.is_open()) {

		error ("could not open file '%s' !",mapfilename[0]);
	}

	if (verbose) {
		Rprintf("Reading map from file '%s' ...\n",mapfilename[0]);
	}

	int mapline=0;
	// read and ignore header line
	if (mapHasHeaderLine) getline(mapfile,data);
	if (strandid==3) {
		while (getline(mapfile,data)) {
			istringstream datas (data);
			mapline++;
			//      if ( datas >> tmp_chrom >> tmp_snpnm >> tmp_phymap >> tmp_strand >> tmp_coding) {
			if ( datas >> tmp_chrom >> tmp_snpnm >> tmp_phymap >> tmp_strand) {
				chrom.push_back(tmp_chrom);
				snpnm.push_back(tmp_snpnm);
				phymap.push_back(tmp_phymap);
				if (tmp_strand=="+") strand.push_back(1);
				else if (tmp_strand=="-") strand.push_back(2);
				else if (tmp_strand=="u") strand.push_back(0);
				else error ("Strand code not recognised at line %i !\n", mapline);
				//	coding.push_back(tmp_coding);
			} else {
				Rprintf("%s\n",tmp_chrom.c_str());
				Rprintf("%s\n",tmp_snpnm.c_str());
				Rprintf("%lu\n",tmp_phymap);
				Rprintf("%s\n",tmp_strand.c_str());
				error ("incomplete map record in line %i:\n%s\n",mapline,data.c_str());
			}
		}
	} else {
		while (getline(mapfile,data)) {
			istringstream datas (data);
			mapline++;
			if ( datas >> tmp_chrom >> tmp_snpnm >> tmp_phymap) {
				chrom.push_back(tmp_chrom);
				snpnm.push_back(tmp_snpnm);
				phymap.push_back(tmp_phymap);
				if (strandid==0) strand.push_back(0);
				else if (strandid==1) strand.push_back(1);
				else if (strandid==2) strand.push_back(2);
			} else {
				error ("incomplete map record in line %i\n",mapline);
			}
		}
	}

	mapfile.close();

	///////////////////////////
	// reading pedigree file //
	///////////////////////////

	int nsnps = chrom.size();

	if (verbose) {
		Rprintf("... done.  Read positions of %i markers from file '%s'\n",nsnps,mapfilename[0]);
	}

	char* chgt[2*MAXIDS];

	ifstream pedfile (pedfilename[0]);
	if( !pedfile.is_open()) {

		error ("could not open file '%s' !",pedfilename[0]);
	}

	if (verbose) {
		Rprintf("Reading genotypes from file '%s' ...\n",pedfilename[0]);
	}

	linecount = 0;
	double lasti = 1;
	if (format==0) {
		while (getline(pedfile,data)) {
			data = replace_mrl(data);
			++linecount;
			if (linecount >= MAXIDS) error ("Number of people is greater than MAXIDS allowed (%i).\nContact GenABEL supporter.\n",MAXIDS);

			istringstream datas (data);


			if (datas >> token >> tmp_iid >> token >> token >> token) {
				//char gdata;
				for (int i=0;i<traits;i++) {datas >> token;}
				iid.push_back(tmp_iid);
				if (!(chgt[2*(linecount-1)] = new char [2*nsnps])) error ("ran out of memory ...\n");
				if (!(chgt[2*(linecount-1)+1] = new char [2*nsnps])) error ("ran out of memory ...\n");
				for (int snp = 0; snp < nsnps; snp++) {
					if (datas >> gdata) {chgt[2*(linecount-1)][snp] = gdata;}
					else {
						error ("too few genotypes for person '%s' ('%s', line %li) !\n",tmp_iid.c_str(),pedfilename[0],linecount);
					}
					if (datas >> gdata) {chgt[2*(linecount-1)+1][snp] = gdata;}
					else {
						error ("too few genotypes for person '%s' ('%s', line %li) !\n",tmp_iid.c_str(),pedfilename[0],linecount);
					}
				}
			} else {
				error("too few records at line %li of pedigree file '%s' !\n",linecount,pedfilename[0]);
			}

			if (verbose && (nsnps*linecount) > ((*bcast)*lasti)) {
				Rprintf("  ... read %li genotypes ...\n",(nsnps*linecount));
				lasti += 1.;
			}

		}
	} else {
		while (getline(pedfile,data)) {
			data = replace_mach(data);
			data = replace_mrl(data);
			++linecount;

			istringstream datas (data);

			//char gdata;

			if (datas >> token >> tmp_iid >> token) {
				iid.push_back(tmp_iid);
				if (!(chgt[2*(linecount-1)] = new char [2*nsnps])) error ("ran out of memory ...\n");
				if (!(chgt[2*(linecount-1)+1] = new char [2*nsnps])) error ("ran out of memory ...\n");
				for (int snp = 0; snp < nsnps; snp++) {
					if (datas >> gdata) {chgt[2*(linecount-1)][snp] = gdata;}
					else {
						error ("too few genotypes for person '%s' ('%s', line %li) !\n",tmp_iid.c_str(),pedfilename[0],linecount);
					}
					if (datas >> gdata) {chgt[2*(linecount-1)+1][snp] = gdata;}
					else {
						error ("too few genotypes for person '%s' ('%s', line %li) !\n",tmp_iid.c_str(),pedfilename[0],linecount);
					}
				}
				//	      	Rprintf("%s ",tmp_iid.c_str());
			} else {
				error("too few records at line %li of pedigree file '%s' !\n",linecount,pedfilename[0]);
			}

			if (verbose && (nsnps*linecount) > ((*bcast)*lasti)) {
				Rprintf("  ... read %li genotypes ...\n",(nsnps*linecount));
				lasti += 1.;
			}

		}
	}
	pedfile.close();

	if (verbose) {
		Rprintf("...done.  Read information for %li people from file '%s'\n",iid.size(),pedfilename[0]);
	}

	if (verbose) {
		Rprintf("Analysing marker information ...\n");
	}

	int nids = iid.size();
	int nbytes = (int)ceil((double)nids/4.);

	int idx;
//	char gdata;
	int* gnum = new int [nids*2];

	int byte;
	int ind;
	int offset[4] = {6,4,2,0};

	lasti = 1.;

	for (int snp = 0; snp < nsnps; snp++) {

		char allele1 = 0;
		char allele2 = 0;


		unsigned long int ca1 = 0;
		unsigned long int ca2 = 0;
		for (idx = 0; idx < 2*nids; idx++) {

			gdata = chgt[idx][snp];
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
					ca1++;
				} else if (allele2 == 0) {
					allele2 = gdata;
					gnum[idx] = 3;
					ca2++;
				} else {
					error ("illegal genotype (three alleles) snp '%s' file '%s' line %li !",
							snpnm[snp].c_str(),pedfilename[0],linecount);
				}
				//	      }
			}

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
					error ("illegal genotype (half missing) snp '%s' file '%s' line %li !",
							snpnm[snp].c_str(),pedfilename[0],linecount);
				}
				idx += 2;
				if (idx >= 2*nids) break;
			}
		}

		gtype.push_back(tmp_gtype);


		if (verbose && ((1+snp)*nids) > ((*bcast)*lasti)) {
			Rprintf("  ... analysed %li genotypes ...\n",((1+snp)*nids));
			lasti += 1.;
		}

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

	for (int i=0;i<nsnps;i++) {
		outfile.width(2);
		outfile.fill('0');
		outfile << intcoding[i] << " ";
	}
	outfile << endl;

	for (int i=0;i<nsnps;i++) {
		outfile.width(2);
		outfile.fill('0');
		outfile << strand[i] << " ";
	}
	outfile << endl;

	for (int i=0;i<nsnps;i++) {
		tmp_gtype = gtype[i];

		for (byte = 0; byte < nbytes; ++byte) {
			outfile.width(2);
			outfile.fill('0');
			outfile << (unsigned int)tmp_gtype[byte];
			outfile << " ";
		}
		outfile << endl;

	}

		delete [] gnum;
		delete [] tmp_gtype;

		if (verbose) {
		Rprintf("... done.\n");
	}

}
}
