#include <cstdlib>  
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std; 

//STL
#include <list>
#include <vector>
#include <map>
#include <iterator>

#include <R.h>




extern "C" {
  void convert_snp_tped (char** tpedfilename, char** tfamfilename, char** outfilename, int* bcast) {

    int verbose = *bcast ? 1 : 0;

    long int linecount;
    string data;
    string token;
    
    list<string> iid; string tmp_iid;
    list<int> chrom; int tmp_chrom;
    list<string> snpnm; string tmp_snpnm;
    list<double> genmap; double tmp_genmap;
    list<unsigned long> phymap; unsigned long tmp_phymap;
    list<unsigned char*> gtype; unsigned char* tmp_gtype;

    ///////////////////////
    // read the tfamfile //
    ///////////////////////

    ifstream tfamfile (tfamfilename[0]);
    if (tfamfile == NULL) {
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
    if (tpedfile == NULL) {
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

	char allele1 = 0;
	char allele2 = 0;

	for (idx = 0; idx < 2*nids; ++idx) {
	  if (datas >> gdata) {

	    if (gdata == allele1) {
	      gnum[idx] = 1;
	    } else if (gdata == allele2) {
	      gnum[idx] = 3;
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
	      tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)1 << offset[ind]);
	      break;
	    case 4:
	      tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)2 << offset[ind]);
	      break;
	    case 6:
	      tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)3 << offset[ind]);
	      break;
	    case 0:
	      // tmp_gtype[byte] = tmp_gtype[byte] | (0 << offset[ind]); // this does nothing
	      break;
	    default:
	      error ("illegal genotype (half missing) snp '%s' file '%s' line %li !",
		     tmp_snpnm.c_str(),tpedfilename[0],linecount);
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
    if (outfile == NULL) {
      error ("could not open file '%s' !",outfilename[0]);
    }

    if (verbose) {
      Rprintf("Writing to file '%s' ...\n",outfilename[0]);
    }    
    
    copy(iid.begin(), iid.end(), ostream_iterator<string>(outfile, " "));
    outfile << endl;

    copy(snpnm.begin(), snpnm.end(), ostream_iterator<string>(outfile, " "));
    outfile << endl;

    copy(chrom.begin(), chrom.end(), ostream_iterator<int>(outfile, " "));
    outfile << endl;

    copy(phymap.begin(), phymap.end(), ostream_iterator<unsigned long>(outfile, " "));
    outfile << endl;

    outfile.flags(hex);

    do {
      tmp_gtype = gtype.front();
      
      for (byte = 0; byte < nbytes; ++byte) {
	outfile << (unsigned int)tmp_gtype[byte];
	outfile << " ";
      }
      outfile << endl;

      delete [] tmp_gtype;

      gtype.pop_front();
    } while (!gtype.empty());
   
    if (verbose) {
      Rprintf("... done.\n");
    }    
 
  }
}



