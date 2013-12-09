#include <new>
#include "gwaa_cpp.h"



#ifdef __cplusplus
extern "C" {
#endif


	/*
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000

//
// Written by Jan Wigginton
	 */
	double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
	{
		if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
		{
			//printf("FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a"
			//		" negative count", obs_hets, obs_hom1, obs_hom2);
			//exit(EXIT_FAILURE);
			error("FATAL ERROR - SNP-HWE: Current genotype configuration includes a negative count");
		}

		int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
		int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

		int rare_copies = 2 * obs_homr + obs_hets;
		int genotypes   = obs_hets + obs_homc + obs_homr;

		double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
		if (het_probs == NULL)
		{
			//printf("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities" );
			//exit(EXIT_FAILURE);
			error("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities");
		}

		int i;
		for (i = 0; i <= rare_copies; i++)
			het_probs[i] = 0.0;

		// start at midpoint
		int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

		// check to ensure that midpoint and rare alleles have same parity
		if ((rare_copies & 1) ^ (mid & 1))
			mid++;

		int curr_hets = mid;
		int curr_homr = (rare_copies - mid) / 2;
		int curr_homc = genotypes - curr_hets - curr_homr;

		het_probs[mid] = 1.0;
		double sum = het_probs[mid];
		for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
		{
			het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
                            				   / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
			sum += het_probs[curr_hets - 2];

			// 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
			curr_homr++;
			curr_homc++;
		}

		curr_hets = mid;
		curr_homr = (rare_copies - mid) / 2;
		curr_homc = genotypes - curr_hets - curr_homr;
		for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
		{
			het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
					/((curr_hets + 2.0) * (curr_hets + 1.0));
			sum += het_probs[curr_hets + 2];

			// add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
			curr_homr--;
			curr_homc--;
		}

		for (i = 0; i <= rare_copies; i++)
			het_probs[i] /= sum;

		// alternate p-value calculation for p_hi/p_lo
   //double p_hi = het_probs[obs_hets];
   //for (i = obs_hets + 1; i <= rare_copies; i++)
   //  p_hi += het_probs[i];

   //double p_lo = het_probs[obs_hets];
   //for (i = obs_hets - 1; i >= 0; i--)
   //   p_lo += het_probs[i];


   //double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;


		double p_hwe = 0.0;
		//  p-value calculation for p_hwe
		for (i = 0; i <= rare_copies; i++)
		{
			if (het_probs[i] > het_probs[obs_hets])
				continue;
			p_hwe += het_probs[i];
		}

		p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

		free(het_probs);

		return p_hwe;
	}


	void snp_summary_exhweWrapper(double *indata, unsigned long int indataHeight,
			unsigned long int indataWidth,	double *outdata,
			unsigned long int &outdataNcol, unsigned long int &outdataNrow,
			unsigned int narg, SEXP *argList)
	{
		//Rprintf("indataHeight=%d\n",indataHeight);
		//Rprintf("indataWidth=%d\n",indataWidth);
		//unsigned int gt[indataHeight];
		unsigned int * gt = new (std::nothrow) unsigned int [indataHeight];
		if (gt == NULL) {
			Rprintf("cannot get RAM for gt\n");
			return;
		}
		if (indata) {
			unsigned int i;
			for(i=0;i<indataHeight*indataWidth;i++) {
				if (ISNAN(indata[i])) gt[i] = 0;
				else gt[i] = 1 + (unsigned int) indata[i];
				//Rprintf(" %f %d;",indata[i],gt[i]);
			}
			//Rprintf("\n");

			unsigned int nids = indataHeight*indataWidth;

			//Rprintf("nids=%d\n",nids);
			snp_summary_exhwe_Processor(gt, nids, outdata);
			//for (int i=0;i<9;i++) Rprintf(" %f",outdata[i]); Rprintf("\n");
		} else {
			outdataNcol = 9;
			outdataNrow = 1;
		}
		delete [] gt;
	}

	void snp_summary_exhwe_Processor(unsigned int *gt, unsigned int nids, double *out) {
		unsigned int i; //,j,idx;
		//unsigned int nids = (*Nids);
		//char str;
		unsigned int count[3];
		double meaids,p;
		count[0]=count[1]=count[2]=0;
		p = 0.;
		for (i=0;i<9;i++) out[i] = 0.;
		for (i=0;i<nids;i++)
			if (gt[i]) {
				count[gt[i]-1]++;
				p+=(gt[i]-1);
			}

		meaids = 1.*(count[0]+count[1]+count[2]);
		out[0] = meaids;
		out[1] = meaids/nids;
		if (meaids>0)
			out[2] = p/(2.*meaids);
		else
			out[2] = 0.0;
		out[3] = count[0];
		out[4] = count[1];
		out[5] = count[2];
		if (meaids>0) {
			double qmax, maf, pmax, loglik0, loglik1, chi2lrt, fmax;
			out[6] = SNPHWE(count[1],count[0],count[2]);
			pmax = out[2];
			qmax = 1.-pmax;
			maf = qmax; if (pmax<qmax) maf = pmax;
			if (maf>1.e-16) {
				fmax = (4.*count[0]*count[2] - 1.*count[1]*count[1])/((2.*count[0]+1.*count[1])*(2.*count[2]+1.*count[1]));
				loglik0 = 0.;
				if (count[0]) loglik0 += 2.*count[0]*log(qmax);
				if (count[1]) loglik0 += 1.*count[1]*log(2.*qmax*pmax);
				if (count[2]) loglik0 += 2.*count[2]*log(pmax);
				loglik1 = 0.;
				if (count[0]) loglik1 += 1.*count[0]*log(qmax*qmax+qmax*pmax*fmax);
				if (count[1]) loglik1 += 1.*count[1]*log(2.*qmax*pmax*(1.-fmax));
				if (count[2]) loglik1 += 1.*count[2]*log(pmax*pmax+qmax*pmax*fmax);
				chi2lrt = 2*(loglik1-loglik0);
				out[7] = fmax;
				out[8] = chi2lrt;
			} else {
				out[7] = 0.;//maf;
				out[8] = 0.;
			}
		} else {
			out[6] = 1.0;
		}

	}


#ifdef __cplusplus
}
#endif

