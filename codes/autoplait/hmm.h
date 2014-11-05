/* 
 * hmm.h --- written by Yasuko Matsubara 2012.2.29
 */
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>

#define ZERO    0.000001 //0.00000001 
#define ONE     0.999999 //0.99999999 
#define INF     1000000000
#define MAXITER 10 //20 
#define BUFSIZE 254 

#define VARMAX  INF 
#define VARMIN  ZERO 

#define MAXVAL 1.0 


#define PI M_PI //3.141592653589
#define E  M_E //2.71828

 
typedef struct {
	int n;		 // # of sequences
	int k;		 // # of states
	double	*pi;	 // initial vector pi[k] 
	double	**A;	 // transition matrix: A[k][k]
	double pi_denom; // pi_denom
	double *A_denom; // A_denom[k]  ... for incremental estimation

	// numerical
	int d;		 // # dimension
	double **mean;	 // GMM mean: 	 mean[k][d]
	double **var;	 // GMM variance: var[k][d]
	double **sum_w;  // total weighted count count[k][d]
	double **M2;	 // M2 for increlemtal var computation : M2[k][d] 
} HMM;


typedef struct {
	// label
	int   id;
	char  tag[BUFSIZE]; 
	int   parent;
	int   st;
	// data
	float **O;	//O[m][d]
	int   m;
	// pattern
	int pid;
}Input;

    
typedef struct {
        double  **alpha;   //[mx][k]
        double  **beta;    //[mx][k]
        double  ***gamma;  //[n][mx][k]
        double  **gamma_space;  //[n*mx][k]
	double  ****xi;	   //[n][mx][k][k]
	double  *scale;    //[mx]	
	int	**idx;	   //[n][mx]; //for k-means clustering
	HMM 	chmm;
} BAUM;

typedef struct {
        double  **delta;   	//[m][k]
        int	**psi;    	//[m][k]
        int	*q;		//[m]
	double  *piL;		//[k]
	double  **AL;		//[k][k]
	double  **biot;		//[k][m]
} VITERBI;



/// hmmutils.c
void LoadHMM(char *fn, HMM *phmm);
void ReadHMM(FILE *fp, HMM *phmm);
void PrintHMM(FILE *fp, HMM *phmm);
void InitHMM(HMM *phmm, int K, int D);
void ResetHMM(HMM *phmm, int K, int D);
void FreeHMM(HMM *phmm);
void SubstHMM(HMM *phmm1, HMM *phmm2);
//
void GenSequenceArray(HMM *phmm, int m, float **O, int *q);
int  GenInitalState(HMM *phmm);
int  GenNextState(HMM *phmm, int q_t);
int  GenSymbol(HMM *phmm, int q_t);
void GenGaussian(HMM *phmm, int q_t, float *O);
//
double pdfL(HMM *phmm, int kid, float *O);
double pdf(HMM *phmm, int kid, float *O);
double GaussianPDF(double mean, double var, double x);


/// forbackward.c 
double Forward(HMM *phmm, float **O, int m, double **alpha, double *scale);
double Backward(HMM *phmm, float  **O, int m, double **beta, double *scale);

/// baum.c
double BaumWelch(HMM *phmm, int n, Input *xlst, BAUM *baum, int wantKMEANS);
void AllocBaum(BAUM *baum, int n, int m, int k, int d);
void FreeBaum(BAUM *baum, int n, int m, int k);

/// viterbi.c
void   AllocViterbi(VITERBI *viter, int n, int m, int k);
double Viterbi(HMM *phmm, int m, float **O, VITERBI *vit);
double ViterbiL(HMM *phmm, int m, float **O, VITERBI *vit);
void FreeViterbi(VITERBI *vit, int n, int m, int k);

/// kmeans.c
void Kmeans();


#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
 

/// dynamic.c
void Split(Input *X, int st, int len, Input *x);


