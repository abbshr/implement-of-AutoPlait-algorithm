/* 
 * baum.c --- written by Yasuko Matsubara 2012.2.29
 */

#include <stdio.h> 
#include <math.h>
#include "nrutil.h"
#include "hmm.h"
#include "tool.h"


#define DELTA 1 
#define DBG 0

/// constant pi values (uniform distrib.)
#define PIC 0  // default: 0 (NO)
void _computeParams();
void _computeGamma();
void _computeXi();

/// both for batch & incremental estimation
double BaumWelch(HMM *phmm, int n, Input *xlst, BAUM *baum, int wantKMEANS){
	int	i, j, k, t, l, r;

	l = 0;
	double **alpha  = baum->alpha;
	double **beta   = baum->beta;
	double ***gamma = baum->gamma;
	double ****xi   = baum->xi;
	double *scale   = baum->scale;
	double delta;
	double Lpreb, Lsum, Lf, Lb;
	// for kmeans
	int **idx = baum->idx;

	if(n==0)
	  error("estimation error! n==0", "BaumWelch");
	// if k==1, nothing to do
	if(k==1){
		Kmeans(phmm, n, xlst, phmm->d, phmm->k, idx);
		return -1;
	}
	if(wantKMEANS){
		// if initial stage and want k-means
		if(phmm->n==0){
			#if(DBG)
			fprintf(stderr,"running k-means...\n");
			#endif
			ResetHMM(phmm, phmm->k, phmm->d);
			Kmeans(phmm, n, xlst, phmm->d, phmm->k, idx);
		}//if
	}//wantKMEANS

	#if(DBG)
	fprintf(stderr,"running baumwelch...\n");
	if(phmm->n > 0)
		fprintf(stderr, "BaumWelch: ++incremental inference\n");
	if(phmm->n > 0 && n < 0)
		fprintf(stderr, "BaumWelch: --decremental inference\n");
	#endif
	n = (int)fabs((float)n);

	int prev_n = phmm->n;
	phmm->n  = n;	 
	SubstHMM(phmm, &baum->chmm);

	// for each sequence
        Lsum = 0.0;
	for(r=0;r<phmm->n;r++){
		float **O = xlst[r].O;
		int m = xlst[r].m;
		Lf = Forward( phmm, O, m, alpha, scale); 
		Lb = Backward(phmm, O, m, beta,  scale); 
		_computeGamma(phmm, alpha, beta, gamma[r], m);
		_computeXi(phmm, O, alpha, beta, xi[r], m);
		Lsum += Lf;
	}//endfor
	// log log P(O |model)
	Lpreb = Lsum; 

	do  {	
		#if(DBG)
		fprintf(stderr, "%d ", l);
		#endif
		///
		/// M-STEP
		///

		/// compute denoms & numerators
		SubstHMM(&baum->chmm, phmm);
		_computeParams(phmm, gamma, xi, xlst);

		///
		/// E-STEP
		///
		Lsum = 0.0;
		for(r=0;r<phmm->n;r++){
			float **O = xlst[r].O;
			int m = xlst[r].m;
			Lf = Forward( phmm, O, m, alpha, scale); 
			Lb = Backward(phmm, O, m, beta,  scale); 
			_computeGamma(phmm, alpha, beta, gamma[r], m);
			_computeXi(phmm, O, alpha, beta, xi[r], m);
			Lsum += Lf;
		}//r
		delta = Lpreb - Lsum;
		Lpreb = Lsum;
		l++;
	}
	//while (l < MAXITER); 
	while (fabs(delta) > DELTA && l < MAXITER); 

	// for incremental fitting
	_computeParams(phmm, gamma, xi, xlst);

	/// avoid error
	if(isnan(Lsum)){
		fprintf(stderr, "baumWelch: isnan, resetHMM...\n");
		ResetHMM(phmm, phmm->k, phmm->d);
	}

	phmm->n = phmm->n + prev_n;
	#if(DBG)
	fprintf(stderr,"baumwelch: Lh=%e\n", Lsum);
	#endif
	return Lsum;

}

void _computeParams(HMM *phmm, double ***gamma, double ****xi, Input *xlst){
	int	r, i, j, k, t;
	/// initial prob vector pi
	// (a) recover previous pi*N 
	for(i=0;i<phmm->k;i++)
		phmm->pi[i] *= phmm->pi_denom;
	// (b) add new gamma (f=1: increment, f=-1: decrement)
	for(i=0;i<phmm->k;i++)
		for(r=0;r<phmm->n;r++)
			phmm->pi[i] += (ZERO+gamma[r][0][i]);
	for(i=0;i<phmm->k;i++)
		if(phmm->pi[i] < 0)
			phmm->pi[i] = 0.0; 
	// (c) normalize
	phmm->pi_denom = 0.0;
	for(i=0;i<phmm->k;i++)
		phmm->pi_denom += phmm->pi[i]; 
	#if(PIC)
	for(i=0;i<phmm->k;i++)
		  phmm->pi[i] = 1.0/(double)phmm->k;
	#else
	for(i=0;i<phmm->k;i++)
		  phmm->pi[i] = phmm->pi[i]/phmm->pi_denom;
	#endif

	/// transition matrix A
	// (a) recover previous A*N 
	for(i=0;i<phmm->k;i++)
		for(j=0;j<phmm->k;j++)
			phmm->A[i][j] *= phmm->A_denom[i];
	// (b) add new xi 
	for(i=0;i<phmm->k;i++)
		for(j=0;j<phmm->k;j++)
			for(r=0;r<phmm->n;r++)
				for(t=0;t<xlst[r].m-1;t++) 
					phmm->A[i][j] += ZERO+xi[r][t][i][j];
	for(i=0;i<phmm->k;i++)
		for(j=0;j<phmm->k;j++)
			if(phmm->A[i][j]<0)
				phmm->A[i][j]=0.0;
	// (c) normalize
	for(i=0;i<phmm->k;i++){
		phmm->A_denom[i] = 0.0;
		for(j=0;j<phmm->k;j++)
			phmm->A_denom[i]+=phmm->A[i][j];
	}//i
	for(i=0;i<phmm->k;i++)
		for(j=0;j<phmm->k;j++)
			phmm->A[i][j] = phmm->A[i][j]/phmm->A_denom[i];

	/// weighted incremental computation
	for(i=0;i<phmm->k;i++){
		for(k=0;k<phmm->d;k++){
			// if initial stage
			if(phmm->sum_w[i][k]  == 0){
				phmm->mean[i][k]   = 0.0;
				phmm->M2[i][k]     = 0.0;
			}//else incremental computation

			for(r=0;r<phmm->n;r++){
				for(t=0;t<xlst[r].m;t++){
					double x = xlst[r].O[t][k];
					double w = gamma[r][t][i];
					double tmp = w + phmm->sum_w[i][k] + ZERO;
					double delta = x - phmm->mean[i][k];
					double R = (delta * w) / tmp;
					phmm->mean[i][k] += R;
					phmm->M2[i][k]   += phmm->sum_w[i][k]*delta*R;
					phmm->sum_w[i][k] = tmp;
				}//t
			}//r
			phmm->var[i][k] = phmm->M2[i][k]/phmm->sum_w[i][k];
			//int N = phmm->m * phmm->n;
			//phmm->var[i][k] *= N/(N-1);
			//double varmn = VARMIN/(double)(phmm->k);
			double varmn = VARMIN;
			// variance setting
			if(phmm->var[i][k] > VARMAX)
				phmm->var[i][k] = VARMAX;
			if(phmm->var[i][k] < varmn)
				phmm->var[i][k] = varmn;


		}//k
	}//i		


}
void _computeGamma(HMM *phmm, double **alpha, double **beta, double **gamma, int m){
	int 	i, j, t;
	for(t=0;t<m;t++){
		double sum= 0.0;
		for(j=0;j<phmm->k;j++){
			gamma[t][j] = alpha[t][j]*beta[t][j];
			sum += gamma[t][j];
		}//j
		for(i=0;i<phmm->k;i++) 
			gamma[t][i] = gamma[t][i]/sum;
	}//t
}
void _computeXi(HMM *phmm, float **O, double **alpha, double **beta, double ***xi, int m){
	int i, j, t;
	for(t=0;t<m-1;t++) {
		double sum = 0.0;	
		for(i=0;i<phmm->k;i++) 
			for(j=0;j<phmm->k;j++){
				xi[t][i][j] = (alpha[t][i]*beta[t+1][j])
					*(phmm->A[i][j] * pdf(phmm, j, O[t+1]));
				sum += xi[t][i][j];
			}//j
		for(i=0;i<phmm->k;i++) 
			for(j=0;j<phmm->k;j++) 
				xi[t][i][j] /= sum;
	}//t
}


void AllocBaum(BAUM *baum, int n, int m, int k, int d){
	int t,r;
	#if(DBG)
	fprintf(stderr, "[AllocBaum] constant pi PIC = %d\n", PIC);	
	#endif
        baum->alpha = dmatrix(0, m, 0, k);
        baum->beta  = dmatrix(0, m, 0, k);
        baum->scale = dvector(0, m);
        baum->idx   = imatrix(0, n, 0, m);
	//gamma
	baum->gamma = (double ***) malloc(n*sizeof(double **));
	for(r=0;r<n;r++)
		baum->gamma[r] = dmatrix(0, m, 0, k);
	
	//xi
	baum->xi = (double ****) malloc(sizeof(double ***)*n);
	for(r=0;r<n;r++){
		baum->xi[r] = (double ***) malloc(sizeof(double **)*m);
		for(t=0;t<m;t++)
			baum->xi[r][t] = dmatrix(0, k, 0, k);
	}//r
	//for incremental EM
	InitHMM(&baum->chmm, k, d);

}
void FreeXi(double *** xi, int m, int k){
	int t;
	for (t = 0; t < m; t++)
		free_dmatrix(xi[t], 0, k, 0, k);
	free(xi);
}
void FreeGamma(double ** gamma, int m, int k){
	free_dmatrix(gamma, 0, m, 0, k);
}
void FreeBaum(BAUM *baum, int n, int m, int k){
	int r;
	for(r=0;r<n;r++){
		FreeXi(baum->xi[r], m, k);
		FreeGamma(baum->gamma[r], m, k);
	}//r
	free_dmatrix(baum->alpha, 0, m, 0, k);
	free_dmatrix(baum->beta, 0, m, 0, k);
	free_dvector(baum->scale, 0, m);
}

