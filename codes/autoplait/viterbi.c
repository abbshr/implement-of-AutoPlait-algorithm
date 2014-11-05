/* 
 * viterbi.c --- written by Yasuko Matsubara 2012.3.1
 */

#include <math.h>
#include "hmm.h"
#include "tool.h"
#include "nrutil.h"

#define VITHUGE  100000000000.0
#define DBG 0 //1

void AllocViterbi(VITERBI *vit, int n, int m, int k){
	vit->delta = dmatrix(0, m, 0, k);
	vit->psi   = imatrix(0, m, 0, k);
	vit->q     = ivector(0, m);
	vit->piL   = dvector(0, k);
	vit->AL    = dmatrix(0, k, 0, k);
	vit->biot  = dmatrix(0, k, 0, m);
}
void FreeViterbi(VITERBI *vit, int n, int m, int k){
	free_dmatrix(vit->delta, 0, m, 0, k);
	free_imatrix(vit->psi, 0, m, 0, k);
	free_ivector(vit->q, 0, m);
	free_dvector(vit->piL, 0, k);
	free_dmatrix(vit->AL, 0, k, 0, k);
	free_dmatrix(vit->biot, 0, k, 0, m);
}

double Viterbi(HMM *phmm, int m, float **O, VITERBI *vit){
	int 	i, j, t;
	int	maxvalind;
	double	maxval, val;
	double **delta = vit->delta;
	int    **psi   = vit->psi;
	int    *q      = vit->q;

	// compute delta (t=0)
	for(i=0;i<phmm->k;i++){
		delta[0][i] = phmm->pi[i] * pdf(phmm, i, O[0]); 
		psi[0][i] = 0;
	}//i	
	// compute delta (t>0)
	for(t=1;t<m;t++){
		for(j=0;j<phmm->k;j++){
			maxval = 0.0;
			maxvalind = 0;	
			for(i=0;i<phmm->k;i++){
				val = delta[t-1][i]*(phmm->A[i][j]);
				if(val>maxval){
					maxval = val;	
					maxvalind = i;	
				}
			}//i
			delta[t][j] = maxval*pdf(phmm, i, O[t]); 
			psi[t][j] = maxvalind; 
		}//j
	}//t
	// final likelihood
	double Lh = 0.0;
	q[m-1] = -1;
	for(i=0;i<phmm->k;i++)
                if(delta[m-1][i]> Lh) {
			Lh = delta[m-1][i];	
			q[m-1] = i;
		}
	// avoid error
	if(q[m-1]==-1)
		error("error", "cannot compute viterbi path");
	//check path
	for(t=m-2;t>=0;t--)
		q[t] = psi[t+1][q[t+1]];

	return Lh;
}

// take logarithm
void TakelogHMM(HMM *phmm, VITERBI *vit, int m, float **O){
	int i,j,t;
	for(i=0;i<phmm->k;i++) 
		vit->piL[i] = log(phmm->pi[i]+ZERO);
	for(i=0;i<phmm->k;i++) 
		for(j=0;j<phmm->k;j++)
			vit->AL[i][j]   = log(phmm->A[i][j]+ZERO);
	for(i=0;i<phmm->k;i++) 
		for(t=0;t<m;t++)
			vit->biot[i][t] = pdfL(phmm, i, O[t]); 
 
}
double ViterbiLList(HMM *phmm, int m, float **O, VITERBI *vit, double *delta){
        int     i, j, t;
        double  maxval, val;

	TakelogHMM(phmm, vit, m, O);

	double deltax = -VITHUGE;

	// compute delta (t==0) 
	deltax = -VITHUGE;
        for(i=0;i<phmm->k;i++){
                vit->delta[0][i] = vit->piL[i] + vit->biot[i][0];
		if(deltax < vit->delta[0][i])
			deltax = vit->delta[0][i];
        }//i
#if(DBG)
 	fprintf(stdout, "%f \n", deltax);
#endif
	delta[0] = deltax;

	// compute delta (t>0) 
        for(t=1;t<m;t++){
		deltax = -VITHUGE;
                for(j=0;j<phmm->k;j++){
                        maxval = -VITHUGE;
                        for(i=0;i<phmm->k;i++){
                                val = vit->delta[t-1][i] + vit->AL[i][j];
                                if (val > maxval) 
                                        maxval = val;
                        }//i
                        vit->delta[t][j] = maxval + vit->biot[j][t]; 
		if(deltax < vit->delta[t][j])
			deltax = vit->delta[t][j];        
                }//j
#if(DBG)
 		fprintf(stdout, "%f \n", deltax);
#endif	
		delta[t] = deltax;	
        }//t

	for(t=m-1;t>=1;t--)
		delta[t]-=delta[t-1];
#if(DBG)
	for(t=0;t<m;t++)
 		fprintf(stdout, "%f \n", delta[t]);
#endif	
		
	
	return deltax; 

}
double ViterbiL(HMM *phmm, int m, float **O, VITERBI *vit){
        int     i, j, t;
        int     maxvalind;
        double  maxval, val;
	if(m==0) return 0;
	TakelogHMM(phmm, vit, m, O);

	double **delta = vit->delta;
	int    **psi   = vit->psi;
	int    *q      = vit->q;
	double *piL    = vit->piL;
	double **AL    = vit->AL;
	double **biot  = vit->biot;

	double deltax = -VITHUGE;

	// compute delta (t==0) 
	deltax = -VITHUGE;
        for(i=0;i<phmm->k;i++){
                delta[0][i] = piL[i] + biot[i][0];
                psi[0][i] = 0;
		if(deltax < delta[0][i])
			deltax = delta[0][i];
        }//i
#if(DBG)
 	fprintf(stdout, "%f \n", deltax);
#endif
	// compute delta (t>0) 
        for(t=1;t<m;t++){
		deltax = -VITHUGE;
                for(j=0;j<phmm->k;j++){
                        maxval = -VITHUGE;
                        maxvalind=0;
                        for(i=0;i<phmm->k;i++){
                                val = delta[t-1][i] + (AL[i][j]);
                                if (val > maxval) {
                                        maxval = val;
                                        maxvalind = i;
                                }
                        }//i
                        delta[t][j] = maxval + biot[j][t]; 
                        psi[t][j] = maxvalind;
		if(deltax < delta[t][j])
			deltax = delta[t][j];        
                }//j
#if(DBG)
 		fprintf(stdout, "%f \n", deltax);
#endif		
        }//t
	
 	// final likelihood
	double Lh = -VITHUGE;
        q[m-1] = -1;
        for(i=0;i<phmm->k;i++){
                if (delta[m-1][i] > Lh) {
                        Lh = delta[m-1][i];
                        q[m-1] = i;
		}
        }//i
 	// avoid error
	if(q[m-1]==-1){
		fprintf(stderr, "Lh=%f\n", Lh);
		PrintHMM(stderr, phmm);
		error("error", "cannot compute viterbi path (LOG)");
	}
	// check viterbi path
	for(t=m-2;t>=0;t--)
		q[t] = psi[t+1][q[t+1]];
	return Lh;
}





