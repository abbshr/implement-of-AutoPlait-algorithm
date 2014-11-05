/* 
 * kmeans.c --- written by Yasuko Matsubara 2012.4.30
 */

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include <sys/types.h>
#include <unistd.h> 

#define MAIN 0
#define DBG 0

#include "hmm.h"
#include "tool.h"
#define KMAXITER 100
#define MINVAR VARMIN
#define RANDOM 0
void Kmeans();


#if MAIN
int main (int argc, char **argv){
	HMM  	hmm; 	 // model
	float	***O; //O[n][m][d]
	double **mean;	 // GMM mean: 	 mean[k][d]
	double **var;	 // GMM variance: var[k][d]
	
  	if(argc != 6)
    	error("usage: cmd seqfn, n(#seqs), m(len), k(#state) d(#dimension)", argv[0]);
	setseed();

	char *fn = argv[1];
	int n  = atoi(argv[2]); // # of seqs
	int mx = atoi(argv[3]); // seq length
  	int k  = atoi(argv[4]); // # of states
  	int d  = atoi(argv[5]); // # of dimensions
	
	int *mlst = ivector(0, n);	
	O = (float***)malloc(sizeof(float**)*n);
	LoadDB(fn, n, mx, d, O, mlst);

  	//PrintSequence(stderr, m, d, O[0]);
	InitHMM(&hmm, k, d);
	int **idx; //idx[n][m];
	idx = imatrix(0, n, 0, m);
	Kmeans(&hmm, n, mlst, d, O, k, idx);

#if DBG
	PrintHMM(stderr, &hmm);
#endif

}
#endif
double distance(float *a, double *b, int d){
	double dist = 0;
	int i;
	for(i=0;i<d;i++)
		dist += (a[i]-b[i])*(a[i]-b[i]);

	return dist;
}
void calcMean(double *mean, int n, Input *xlst, int d, int **idx, int label){
	int i,j,r; int cnt = 0;
	// if label==-1, k=1, single core
	// for each dimension (j=0...d-1)
	for(j=0;j<d;j++){
		cnt = 0;
		mean[j] = 0;
		for(r=0;r<n;r++)
			for(i=0;i<xlst[r].m;i++)
				if(idx[r][i] == label || label==-1){
					mean[j] += xlst[r].O[i][j];
					cnt++;
				}//if
		mean[j] /= (cnt+ZERO);
		if(cnt==0)
			mean[j] = getrand();
	}//j
}

void calcVar(double *var, int n, Input *xlst, int d, int **idx, int label){
	int i,j,r;
	// if label==-1, k=1, single core
	// for each dimension (j=0...d-1)
	for(j=0;j<d;j++){
		int cnt = 0;
		double sum = 0;
		double sum_sqr = 0;
		for(r=0;r<n;r++)
			for(i=0;i<xlst[r].m;i++)
				if(idx[r][i] == label || label==-1){
					cnt++;
					sum += xlst[r].O[i][j];
					sum_sqr += xlst[r].O[i][j]*xlst[r].O[i][j];
				}//if	
		double mean = sum / (cnt+ZERO);
		var[j] = ZERO + (sum_sqr-sum*mean)/(cnt-1+ZERO);
		if(cnt<=1 || var[j] < MINVAR)
			var[j] = MINVAR;

	}//j

}

int _randInt(int m){
	int r = (int) (getrand()*m);
	if(m==0)
		return 0;
	else
		return r;
}
void Kmeans(HMM *phmm, int n, Input *xlst, int d, int k, int **idx){
	int i,j,r;
	int l;
	double **mean = phmm->mean;
	double **var  = phmm->var;

	// if k==1, nothing to do
	if(k==1){
		calcMean(mean[0], n, xlst, d, idx, -1);
		calcVar(  var[0], n, xlst, d, idx, -1);
		return;
	}//

	#if(RANDOM)
	// random assign
	for(r=0;r<n;r++)
		for(i=0;i<xlst[r].m;i++)
			idx[r][i] = (int) (getrand()*10)%k;
	#else
	// uniform assign
	int cnt=0;
	for(r=0;r<n;r++)
		for(i=0;i<xlst[r].m;i++){
			idx[r][i] = (int) (cnt)%k;
			cnt++;
		}//i
	#endif

#if DBG
	for(r=0;r<n;r++)
		for(i=0;i<xlst[r].m;i++)
			fprintf(stderr, "%d ", idx[r][i]);
	fprintf(stderr, "\n");
#endif

	double prev = INF;
	for(l=0;l<KMAXITER;l++){
		double total = 0;
		for(r=0;r<n;r++)
			for(i=0;i<xlst[r].m;i++){
				double min = INF;
				for(j=0;j<k;j++){
					double dist = distance(xlst[r].O[i],mean[j], d); 
					if(min > dist){
						idx[r][i] = j;
						min = dist;
					}//if
				}//j
				total += min;
			}//i

///*
		/// --- avoid single core --- ///
		for(j=0;j<k;j++){
			int frag=0;
			for(r=0;r<n;r++)
				for(i=0;i<xlst[r].m;i++)
					if(idx[r][i]==j){frag=1; break;}
			if(frag==0){
				r = _randInt(n);
				int m = _randInt(xlst[r].m);
				idx[r][m] = j;
			}
		}//j
		/// --- avoid single core --- ///
//*/	
		for(j=0;j<k;j++)
			calcMean(mean[j], n, xlst, d, idx, j);

		if(prev == total)
			break;
		prev = total;
	}//l

#if DBG
	for(r=0;r<n;r++)
		for(i=0;i<xlst[r].m;i++)
			fprintf(stderr, "%d ", idx[r][i]);
	fprintf(stderr, "\n");
#endif

	for(j=0;j<k;j++){
		calcMean(mean[j], n, xlst, d, idx, j);
		calcVar(  var[j], n, xlst, d, idx, j);
	}//j

}

