//#=====================================================================================
//#
//#       Filename:  dometa.h
//#
//#    Description:  Function for meta analysis. 
//#
//#        Version:  1.0
//#        Created:  06-July-2009
//#       Revision:  none
//#
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================

#ifndef SMV_DOMETA_C_H
#define SMV_DOMETA_C_H


extern "C" {
void dometa_c(double *beta_set1, double *beta_set2,
              double *sebeta_set1, double *sebeta_set2,
              double *lambda_set1, double *lambda_set2,
              unsigned *num,
              double *mbeta,
              double *mse);


}
#endif

