/**
 *
 * 2007.07.20 by Yurii Aulchenko, EMCR
 * based on code conributed by Toby Johnson (convert_snp_tped.cpp)
 *
 * last modified 2007.07.23 by YA
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

string replace_mrl(string in) {
	int pos = 0;
	while ( (pos=in.find("/")) != -1 ) {
		in.erase(pos,1);
		in.insert(pos," ");
	}
	return in;
}

string replace_mach(string in) {
	int pos = 0;
	if ( (pos=in.find("->")) != -1 ) {
		in.erase(pos,2);
		in.insert(pos," ");
	}
	return in;
}

