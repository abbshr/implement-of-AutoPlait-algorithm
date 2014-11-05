/* 
 * forbackward.c --- written by Yasuko Matsubara 2012.2.29
 */
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include "hmm.h"



double Forward(HMM *phmm, float **O, int m, double **alpha, double *scale){
	int	i,j,t;
	double sum;
	double L;
        for (t = 0; t < m; t++)
          scale[t] = 0.0;
	/// init 
	scale[0] = 0.0;	
	for (i = 0; i < phmm->k; i++) {
		alpha[0][i] = phmm->pi[i]* pdf(phmm, i, O[0]); //(phmm->B[i][O[0]]);
		scale[0] += alpha[0][i];
	}//i
	for (i = 0; i < phmm->k; i++) 
		alpha[0][i] /= scale[0]; 
	/// induction
	for (t = 0; t < m - 1; t++) {
		scale[t+1] = 0.0;
		for (j = 0; j < phmm->k; j++) {
			sum = 0.0;
			for (i = 0; i < phmm->k; i++) 
				sum += alpha[t][i]* (phmm->A[i][j]); 

			alpha[t+1][j] = sum*(pdf(phmm, j, O[t+1])); //phmm->B[j][O[t+1]]);
			scale[t+1] += alpha[t+1][j];
		}//j
		for (j = 0; j < phmm->k; j++) 
			alpha[t+1][j] /= scale[t+1]; 
	}//t
	/// termination
	L = 0.0;
	for (t = 0; t < m; t++)
		L += log(scale[t]);

	return L;
}

double Backward(HMM *phmm, float **O, int m, double **beta, double *scale){
        int     i, j, t; 
	double sum;
	/// init 
        for (i = 0; i < phmm->k; i++)
                beta[m-1][i] = 1.0/scale[m-1]; 
	///indution 
        for (t = m - 2; t >= 0; t--) {
                for (i = 0; i < phmm->k; i++) {
			sum = 0.0;
                        for (j = 0; j < phmm->k; j++)
                        	sum += phmm->A[i][j] * 
					pdf(phmm, j, O[t+1])*beta[t+1][j];
					//(phmm->B[j][O[t+1]])*beta[t+1][j];
                        beta[t][i] = sum/scale[t];
                }//i
        }//t
        double L = 0.0;
        for (t = 0; t < m; t++)
                L += log(scale[t]);
	return L;
}
