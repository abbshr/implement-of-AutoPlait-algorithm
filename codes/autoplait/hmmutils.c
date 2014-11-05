/* 
 * hmmutils.c --- written by Yasuko Matsubara 2012.2.29
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "hmm.h"
#include "tool.h"
#define Nsample 10
#define VARINI 10 //0000

//------------------//
// HMM (read, free, init, copy, subst)
//------------------//
void LoadHMM(char *fn, HMM *phmm){
	FILE *fp = fopen(fn,"r");
  	if (fp == NULL)
    		error("cannot open file", fn);
  	ReadHMM(fp, phmm);
  	fclose(fp);
}
void ReadHMM(FILE *fp, HMM *phmm){
	int i, j, k;

	// read params
	if(fscanf(fp, "k= %d\n", &(phmm->k))==EOF)
	  error("error","cannot read HMM");
	if(fscanf(fp, "d= %d\n", &(phmm->d))==EOF)
	  error("error","cannot read HMM");

	// alloc etc.
	InitHMM(phmm, phmm->k, phmm->d);
	// read pi, A, mean, var
	fscanf(fp, "pi:\n");
	for (i = 0; i < phmm->k; i++) 
		fscanf(fp, "%lf", &(phmm->pi[i])); 
	fscanf(fp,"\n");
	fscanf(fp, "A:\n");
	for(i=0;i<phmm->k;i++){ 
		for(j=0;j<phmm->k;j++)
			fscanf(fp, "%lf", &(phmm->A[i][j])); 
		fscanf(fp,"\n");
	}//i
	fscanf(fp, "mean:\n");
	for(j=0;j<phmm->k;j++){ 
		for(k=0;k<phmm->d;k++) 
			fscanf(fp, "%lf", &(phmm->mean[j][k])); 
		fscanf(fp,"\n");
	}//j
	fscanf(fp,"\n");
	fscanf(fp, "var:\n");
	for(j=0;j<phmm->k;j++){ 
		for(k=0;k<phmm->d;k++) 
			fscanf(fp, "%lf", &(phmm->var[j][k])); 
		fscanf(fp,"\n");
	}//j

}
void FreeHMM(HMM *phmm){
	free_dvector(phmm->pi, 0, phmm->k);
	free_dmatrix(phmm->A,  0, phmm->k, 0, phmm->k);
	free_dvector(phmm->A_denom, 0, phmm->k);
	free_dmatrix(phmm->mean, 0, phmm->k, 0, phmm->d);
	free_dmatrix(phmm->var,  0, phmm->k, 0, phmm->d);
	free_dmatrix(phmm->sum_w,0, phmm->k, 0, phmm->d);
	free_dmatrix(phmm->M2,  0, phmm->k, 0, phmm->d);
}
void ResetHMM(HMM *phmm, int K, int D){
	int i, j, k;
	double sum;
	phmm->n = 0;
       	phmm->d = D;
        phmm->k = K;

	// A matrix
        for(i=0;i<phmm->k;i++){
		sum = 0.0;
                for(j=0;j<phmm->k;j++){
                        phmm->A[i][j] = getrand(); 
			sum += phmm->A[i][j];
		}//j
                for(j=0;j<phmm->k;j++) 
			 phmm->A[i][j] /= sum;
	}//i

	// initial pi
	sum = 0.0;
        for(i=0;i<phmm->k;i++){
                phmm->pi[i] = getrand(); 
		sum += phmm->pi[i];
	}//i
        for(i=0;i<phmm->k;i++) 
		phmm->pi[i] /= sum;

	// gaussian params - mean & variance
        for(j=0;j<phmm->k;j++) {
		sum = 0.0;	
                for(k=0;k<phmm->d;k++){
                        phmm->mean[j][k] = getrand()*MAXVAL;
                        phmm->var[j][k]  = VARINI;
		}//k
	}//j

	//for incremental maintenance
	phmm->pi_denom  = 0;
	for(i=0;i<phmm->k; i++){
		phmm->A_denom[i] = 0;
		for(k=0;k<phmm->d;k++){
			phmm->sum_w[i][k] = 0.0;
			phmm->M2[i][k]    = 0.0;
		}//k
	}//i

}
void InitHMM(HMM *phmm, int K, int D){
	int i, j, k;
	double sum;
	phmm->n = 0;
       	phmm->d = D;
        phmm->k = K;

	// A matrix
        phmm->A = (double **) dmatrix(0, phmm->k, 0, phmm->k);
        for(i=0;i<phmm->k;i++){
		sum = 0.0;
                for(j=0;j<phmm->k;j++){
                        phmm->A[i][j] = getrand()+ZERO; 
			sum += phmm->A[i][j];
		}//j
                for(j=0;j<phmm->k;j++) 
			 phmm->A[i][j] /= sum;
	}//i

	// initial pi
        phmm->pi = (double *) dvector(0, phmm->k);
	sum = 0.0;
        for(i=0;i<phmm->k;i++){
                phmm->pi[i] = getrand()+ZERO; 
		sum += phmm->pi[i];
	}//i
        for(i=0;i<phmm->k;i++) 
		phmm->pi[i] /= sum;

	// gaussian params - mean & variance
        phmm->mean = (double **) dmatrix(0, phmm->k, 0, phmm->d);
        phmm->var  = (double **) dmatrix(0, phmm->k, 0, phmm->d);
        for(j=0;j<phmm->k;j++) {
		sum = 0.0;	
                for(k=0;k<phmm->d;k++){
                        phmm->mean[j][k] = getrand()*MAXVAL;
                        phmm->var[j][k]  = VARINI;
		}//k
	}//j

	//for incremental maintenance
	phmm->pi_denom  = 0;
	phmm->A_denom   = (double *)  dvector(0, phmm->k);
	phmm->sum_w     = (double **) dmatrix(0, phmm->k, 0, phmm->d);
	phmm->M2        = (double **) dmatrix(0, phmm->k, 0, phmm->d);
	for(i=0;i<phmm->k; i++){
		phmm->A_denom[i] = 0;
		for(k=0;k<phmm->d;k++){
			phmm->sum_w[i][k] = 0.0;
			phmm->M2[i][k]    = 0.0;
		}//k
	}//i

}
// phmm1:from, phmm2:to
void SubstHMM(HMM *phmm1, HMM *phmm2){
        int i, j, k;
        phmm2->d = phmm1->d;
        phmm2->k = phmm1->k;
        phmm2->n = phmm1->n;
        for(i=0;i<phmm2->k;i++)
                for(j=0;j<phmm2->k;j++)
                        phmm2->A[i][j] = phmm1->A[i][j];
        for(j=0;j<phmm2->k;j++)
                for(k=0;k<phmm2->d;k++){
                        phmm2->mean[j][k] = phmm1->mean[j][k];
                        phmm2->var[j][k]  = phmm1->var[j][k];
			phmm2->sum_w[j][k]= phmm1->sum_w[j][k];
                        phmm2->M2[j][k]  = phmm1->M2[j][k];
		}//k
        for(i=0;i<phmm2->k;i++)
                phmm2->pi[i] = phmm1->pi[i]; 
	//for incremental maintenance
	phmm2->pi_denom = phmm1->pi_denom;
	for (i = 0; i < phmm2->k; i++){
		phmm2->A_denom[i] = phmm1->A_denom[i];
	}//i
}


void PrintHMM(FILE *fp, HMM *phmm){
        int i, j, k;
	fprintf(fp, "k= %d\n", phmm->k); 
	fprintf(fp, "d= %d\n", phmm->d); 
	fprintf(fp, "pi:\n");
        for (i = 0; i < phmm->k; i++) 
		fprintf(fp, "%f ", phmm->pi[i]);
	fprintf(fp, "\n");
	fprintf(fp, "A:\n");
        for (i = 0; i < phmm->k; i++) {
                for (j = 0; j < phmm->k; j++) 
                        fprintf(fp, "%f ", phmm->A[i][j] );
		fprintf(fp, "\n");
	}//i
	fprintf(fp, "mean:\n");
        for (j = 0; j < phmm->k; j++) {
                for (k = 0; k < phmm->d; k++)
                        fprintf(fp, "%f ", phmm->mean[j][k]);
		fprintf(fp, "\n");
	}//j
	fprintf(fp, "var:\n");
        for (j = 0; j < phmm->k; j++) {
                for (k = 0; k < phmm->d; k++)
                        fprintf(fp, "%f ", phmm->var[j][k]);
		fprintf(fp, "\n");
	}//j

}

//------------------//
// HMM generator
//------------------//
void GenSequenceArray(HMM *phmm, int m, float **O, int *q){
        int t;
        int q_t, o_t;
        q[0] = GenInitalState(phmm);
        GenGaussian(phmm, q[0], O[0]);
        for (t = 1; t < m; t++) {
                q[t] = GenNextState(phmm, q[t-1]);
                GenGaussian(phmm, q[t], O[t]);
        }//t
	//PrintSequence(stderr, m, phmm->d, O);
}
int GenInitalState(HMM *phmm){
        int i;
	double cum = 0.0;
        double rand = getrand();
	for(i=0;i<phmm->k;i++){
		cum += phmm->pi[i];
		if(rand < cum + ZERO)
			return i;
	}//i
	PrintHMM(stderr, phmm);
	fprintf(stderr, "rand=%f, cum=%f\n", rand, cum);
	error("error","cannot assign (pi)");
	return -1;
}
int GenNextState(HMM *phmm, int q_prev){
        int j;
        double cum = 0.0;
        double rand = getrand();
	for(j=0;j<phmm->k;j++){
		cum += phmm->A[q_prev][j];
		if(rand < cum + ZERO) 	
			return j;
	}//j
	error("error","cannot assign (A)");        
	return -1;
}

double Gaussian(double mean, double var){
   double rand;
   double X = 0.0;
   int i;
   /// Central Limit Thorem Method
   for (i=0;i<Nsample;i++){
   	rand = getrand();
	X   += rand;
   }//i
   /// for uniform randoms in [0,1], mu = 0.5 and var = 1/12 
   // adjust X so mu = 0 and var = 1 
   X = (X - Nsample*0.5);	// set mean to 0.0 
   X = X * sqrt(12/Nsample);	// adjust variance to 1
   X = mean + sqrt(var)*X;	// set mean & var
   return X;
}
double GaussianPDF(double mean, double var, double x){
	var = fabs(var); 
	double p = exp(-(x-mean)*(x-mean)/(2*var)) / (sqrt(2*PI*var));
	if(p>=1.0) p=ONE;
	if(p<=0.0) p=ZERO;
	return p; 
}
double pdfL(HMM *phmm, int kid, float *O){
	int i;
	double p=0.0;
	for(i=0;i<phmm->d;i++)
		p += log(ZERO+GaussianPDF(phmm->mean[kid][i], phmm->var[kid][i], (double)O[i]));
	if(p<log(ZERO)) p=log(ZERO);
//fprintf(stderr, "%f\n", p);
	return p;
}
double pdf(HMM *phmm, int kid, float *O){
	int i;
	double p = exp( pdfL(phmm, kid, O) );
	//fprintf(stderr, "pl=%f , p=%f\n", pl, p);
	return p; 
}
void GenGaussian(HMM *phmm, int q_t, float *O){
	int i;
	for(i=0;i<phmm->d;i++)
		O[i] = (float) Gaussian(phmm->mean[q_t][i],phmm->var[q_t][i]); 
}

