#include "export_plink.h"


#ifdef __cplusplus
extern "C" {
#endif

std::string* getGenotype(std::string coding, std::string sep)
{
    std::string* Genotype = new (std::nothrow) std::string [4];
    std::string Letter0 = coding.substr(0,1);
    std::string Letter1 = coding.substr(1,1);
    Genotype[0] = "0"+sep+"0";
    Genotype[1] = Letter0+sep+Letter0;
    Genotype[2] = Letter0+sep+Letter1;
    Genotype[3] = Letter1+sep+Letter1;
    return Genotype;
}

SEXP export_plink(SEXP Ids, SEXP Snpdata, SEXP Nsnps, SEXP NidsTotal,
                  SEXP Coding, SEXP From, SEXP To, SEXP Male, SEXP Traits,
                  SEXP Pedfilename, SEXP Plink, SEXP Append)
{
    int from = INTEGER(From)[0];
    int to = INTEGER(To)[0];
    
		
		if(from <1 || from > to) {error("The function SEXP export_plink(SEXP Ids, SEXP Snpdata, SEXP Nsnps, SEXP NidsTotal,... reports: the variable FROM should be >=1 and less then the variable TO.");} //Maksim 

		
		std::vector<unsigned short int> sex;
    sex.clear();
    unsigned short int sx;
    for(int i=(from - 1); i<to; i++) {
        sx = INTEGER(Male)[i];
        if (sx==0) sx=2;
        //Rprintf("%d %d\n",i,sx);
        sex.push_back(sx);
    }
    std::vector<std::string> ids;
    for(unsigned int i=0; i<((unsigned int) length(Ids)); i++)
        ids.push_back(CHAR(STRING_ELT(Ids,i)));

    std::vector<std::string> coding;
    for(unsigned int i=0; i<((unsigned int) length(Coding)); i++)
        coding.push_back(CHAR(STRING_ELT(Coding,i)));

    //Rprintf("0\n");
    unsigned int nsnps = INTEGER(Nsnps)[0];
    int nids = to - from + 1;
    int nidsTotal = INTEGER(NidsTotal)[0];
    int ntraits = INTEGER(Traits)[0];
    bool append = LOGICAL(Append)[0];
    bool plink = LOGICAL(Plink)[0];
    std::string filename = CHAR(STRING_ELT(Pedfilename,0));
    std::ofstream fileWoA;
    int ieq1 = 1;
    char * snpdata = (char *) RAW(Snpdata);

    //	int gtint[nidsTotal];
    int *gtint = new (std::nothrow) int[nidsTotal];

    //Rprintf("nsnps=%d\n",nsnps);
    //Rprintf("nids=%d\n",nids);
    //Rprintf("to=%d\n", to);
    //Rprintf("from=%d\n", from);

    //char gtMatrix[nids][nsnps];
    char **gtMatrix = new (std::nothrow) char*[nids];
    for (int i=0; i<nids; i++) {
        gtMatrix[i] = new (std::nothrow) char[nsnps];
    }

    //Rprintf("1\n");
    std::string* Genotype;
    std::string sep="/";
    int nbytes;

    //Rprintf("nsnps=%d\n",nsnps);
    //Rprintf("nids=%d\n",nids);

    if ((nids % 4) == 0) {
        nbytes = nidsTotal/4;
    }
    else {
        nbytes = ceil(1.*nidsTotal/4.);
    }

    if (plink) sep=" ";

    if (append)
        fileWoA.open(filename.c_str(),std::fstream::app);
    else
        fileWoA.open(filename.c_str(),std::fstream::trunc);

    //Rprintf("A\n");
    for (unsigned int csnp=0; csnp<nsnps; csnp++) {
        // collect SNP data
        get_snps_many(snpdata+nbytes*csnp, &nidsTotal, &ieq1, gtint);
        for (int iii=from-1; iii<to; iii++) {
            //Rprintf(" %d",gtint[iii]);
            gtMatrix[iii-from+1][csnp] = gtint[iii];
        }
        //Rprintf("\n");
    }

    //Rprintf("B\n");
    for (int i=0; i<nids; i++) {
        fileWoA << i+from << " " << ids[i] << " 0 0 " << sex[i];

        for (int j=0; j<ntraits; j++) fileWoA << " " << 0;
        // unwrap genotypes
        for (unsigned int csnp=0; csnp<nsnps; csnp++) {
            Genotype = getGenotype(coding[csnp], sep);
            // figure out the coding
            fileWoA << " " << Genotype[gtMatrix[i][csnp]];
            delete [] Genotype;
        }
        // end unwrap
        fileWoA << "\n";
    }
    //Rprintf("C\n");
    fileWoA.close();

    //Rprintf("oooo!" );
    //for (int i=0; i<10; i++) Rprintf("%d ",sex[i]);
    //Rprintf("oooo!\n" );

    sex.clear();

    for(int i=0; i<nids; i++) {
        delete [] gtMatrix[i];
    }
    delete [] gtMatrix;

    delete [] gtint;

    return R_NilValue;
}


SEXP export_plink_tped(SEXP Snpnames, SEXP Chromosomes, SEXP Map,
                       SEXP Snpdata, SEXP Nsnps, SEXP Nids, SEXP Coding,
                       SEXP Pedfilename, SEXP ExportNumeric)
{
    std::vector<std::string> snpName;
    for(unsigned int i=0; i<((unsigned int) length(Snpnames)); i++)
        snpName.push_back(CHAR(STRING_ELT(Snpnames,i)));

    std::vector<std::string> coding;
    for(unsigned int i=0; i<((unsigned int) length(Coding)); i++)
        coding.push_back(CHAR(STRING_ELT(Coding,i)));

    std::vector<std::string> chromosome;
    for(unsigned int i=0; i<((unsigned int) length(Chromosomes)); i++)
        chromosome.push_back(CHAR(STRING_ELT(Chromosomes,i)));

    std::vector<double> position;
    for(unsigned int i=0; i<((unsigned int) length(Map)); i++)
        position.push_back(REAL(Map)[i]);

    //Rprintf("0\n");
    unsigned int nsnps = INTEGER(Nsnps)[0];
    int nids = INTEGER(Nids)[0];
    bool exportNumeric = LOGICAL(ExportNumeric)[0];
    std::string filename = CHAR(STRING_ELT(Pedfilename,0));
    std::ofstream fileWoA;
    int ieq1 = 1;
    char * snpdata = (char *) RAW(Snpdata);

    //	int gtint[nids];
    int *gtint = new (std::nothrow) int[nids];

    //Rprintf("nsnps=%d\n",nsnps);
    //Rprintf("nids=%d\n",nids);

    //Rprintf("1\n");
    std::string* Genotype;
    std::string sep=" ";
    int nbytes;

    //Rprintf("nsnps=%d\n",nsnps);
    //Rprintf("nids=%d\n",nids);

    if ((nids % 4) == 0) {
        nbytes = nids/4;
    }
    else {
        nbytes = ceil(1.*nids/4.);
    }

    fileWoA.open(filename.c_str(), std::fstream::trunc);

    //Rprintf("A\n");
    for (unsigned int csnp=0; csnp<nsnps; csnp++) {
        // collect SNP data
        get_snps_many(snpdata+nbytes*csnp, &nids, &ieq1, gtint);
        Genotype = getGenotype(coding[csnp], sep);
        fileWoA << chromosome[csnp] << " " << snpName[csnp]
                << " 0 " << (unsigned long int) position[csnp];

        if (!exportNumeric) {
            for (int i=0; i<nids; i++) {
                fileWoA << " " << Genotype[gtint[i]];
            }
        } else {
            for (int i=0; i<nids; i++) {
                if (gtint[i]==0)
                    fileWoA << " NA";
                else
                    fileWoA << " " << (gtint[i]-1);
            }
        }
        fileWoA << "\n";
    delete [] Genotype;
        //Rprintf("\n");
    }
    //Rprintf("B\n");
    /**
       for (int i=0; i<nids; i++) {
       fileWoA << i+from << " " << ids[i] << " 0 0 " << sex[i];
       for (int j=0; j<ntraits; j++) fileWoA << " " << 0;
       // unwrap genotypes
       for (unsigned int csnp=0; csnp<nsnps; csnp++) {
       Genotype = getGenotype(coding[csnp],sep);
       // figure out the coding
       fileWoA << " " << Genotype[gtMatrix[i][csnp]];
       //fileWoA << " x" << Letter0 << Letter1 << Genotype[0] << Genotype[1] << Genotype[2] << Genotype[3];
       }
       // end unwrap
       fileWoA << "\n";
       }
    **/
    //Rprintf("C\n");
    fileWoA.close();

    //Rprintf("oooo!" );
    //for (int i=0; i<10; i++) Rprintf("%d ",sex[i]);
    //Rprintf("oooo!\n" );

    delete [] gtint;

    return R_NilValue;
}

#ifdef __cplusplus
}
#endif
