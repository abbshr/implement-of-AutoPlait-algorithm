/* 
 * cps.c --- written by Yasuko Matsubara 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include "hmm.h"
#include "tool.h"
#include "autoplait.h"

#define DBG    1 
#define S0  1
#define S1 -1


PlaitWS *ws;
double CPSearch();
double _search_aux();
void AllocCPS();
void FreeCPS();
void _printP(double *P, int k){
	int i;
	for(i=0;i<k;i++)
		fprintf(stderr, "%f ", P[i]);
	fprintf(stderr, "\n");
}
void AllocCPS(CPS *cps, int maxk, int maxlen){
	fprintf(stderr, "alloc cut point search...(k:%d,len:%d)\n", maxk, maxlen);	
	cps->Pu = dvector(0, maxk);
	cps->Pv = dvector(0, maxk);
	cps->Pi = dvector(0, maxk);
	cps->Pj = dvector(0, maxk);
	cps->Su = imatrix(0, maxk, 0, maxlen);
	cps->Sv = imatrix(0, maxk, 0, maxlen);
	cps->Si = imatrix(0, maxk, 0, maxlen);
	cps->Sj = imatrix(0, maxk, 0, maxlen);
	cps->nSu = ivector(0, maxk);
	cps->nSv = ivector(0, maxk);
	cps->nSi = ivector(0, maxk);
	cps->nSj = ivector(0, maxk);
}
void FreeCPS(CPS *cps, int maxk, int maxlen){
	free_dvector(cps->Pu, 0, maxk);
	free_dvector(cps->Pv, 0, maxk);
	free_dvector(cps->Pi, 0, maxk);
	free_dvector(cps->Pj, 0, maxk);
	free_imatrix(cps->Su, 0, maxk, 0, maxlen);
	free_imatrix(cps->Sv, 0, maxk, 0, maxlen);
	free_imatrix(cps->Si, 0, maxk, 0, maxlen);
	free_imatrix(cps->Sj, 0, maxk, 0, maxlen);
	free_ivector(cps->nSu, 0, maxk);
	free_ivector(cps->nSv, 0, maxk);
	free_ivector(cps->nSi, 0, maxk);
	free_ivector(cps->nSj, 0, maxk);
}

double CPSearch(SegBox *Sx, SegBox *s0, SegBox *s1, PlaitWS *wsd){
	int i;
	ws = wsd;
	ResetStEd(s0);
	ResetStEd(s1);
	double lh = 0;
	for(i=0;i<Sx->nSb;i++)
		lh += _search_aux(Sx->Sb[i].st, Sx->Sb[i].len, s0, s1);
	return lh;
} 

int _findMax(double *P, int k){
	int i;
	int loc = -1;
	double max = -INF;
	for(i=0;i<k;i++)
		if(max < P[i]){ max = P[i]; loc = i; }
	return loc;	
}
int _copy_path(int *from, int nfrom, int *to){
	int i;
	for(i=0;i<nfrom;i++)
		to[i]=from[i];
	return nfrom;
}
void _reset_npaths(int *nS, int k){
	int i;
	for(i=0;i<k;i++)
		nS[i] = 0;
}
void _print_path(FILE *fp, int *S, int len){
	int i;
	fprintf(fp, "---path---\n");
	for(i=0;i<len;i++)
		fprintf(fp, "%d\n", S[i]);
	fprintf(fp, "----------\n");
}

double _search_aux(int st, int len, SegBox *s0, SegBox *s1){
	int i,j,u,v,t;
	int k0 = s0->model.k;
	int k1 = s1->model.k;
	int currentID;
	///--- setting ---///
	float **O = ws->x[0].O;
	HMM *m0 = &s0->model;
	HMM *m1 = &s1->model;
	double d0 = s0->delta; 
	double d1 = s1->delta; 
	///--- setting ---///
	double *Pu = ws->cps.Pu;
	double *Pv = ws->cps.Pv;
	double *Pi = ws->cps.Pi;
	double *Pj = ws->cps.Pj;
	int **Su = ws->cps.Su;
	int **Sv = ws->cps.Sv;
	int **Si = ws->cps.Si;
	int **Sj = ws->cps.Sj;	
	int *nSu = ws->cps.nSu;
	int *nSv = ws->cps.nSv;
	int *nSi = ws->cps.nSi;
	int *nSj = ws->cps.nSj;
	_reset_npaths(nSu, k0);
	_reset_npaths(nSv, k0);
	_reset_npaths(nSi, k1);
	_reset_npaths(nSj, k1);
	///--- setting ---///
	// for swap
	HMM *hmm;
	int **imat;
	int *ivec;
	double *dvec;
	///--- setting ---///
	if(d0<=0 || d1<=0)
		error("delta == 0","cpsearch");

	// (t == 0);
	t=st;
	for(v=0;v<k0;v++)
		Pv[v] = log(d1) + log(m0->pi[v]+ZERO) + pdfL(m0, v, O[t]);
	for(j=0;j<k1;j++)
		Pj[j] = log(d0) + log(m1->pi[j]+ZERO) + pdfL(m1, j, O[t]);
	//fprintf(stderr, "Pv:"); _printP(Pv, k0);fprintf(stderr, "Pj:");_printP(Pj, k1);
	
	// (t >= 1);
	for(t=st+1;t<st+len;t++){
		///---  Pu(t)  ---///
		// regime-switch
		int maxj = _findMax(Pj, k1);
		for(u=0;u<k0;u++){
			// if switch
			double maxPj = Pj[maxj] + log(d1) + log(m0->pi[u]+ZERO)   + pdfL(m0, u, O[t]);
			//fprintf(stderr, "(t:%d) Pj: %f j:%d\n", t, maxPj, maxj);
			//fprintf(stderr, "%f +%f +%f +%f \n", Pj[maxj], log(d1), log(m0->pi[u]+ZERO), pdfL(m0, u, O[t]));
			// else
			double maxPv = -INF; int maxv;
			for(v=0;v<k0;v++){
				double val = log(1.0-d0) + Pv[v] + log(m0->A[v][u]+ZERO)  + pdfL(m0, u, O[t]);
				//fprintf(stderr, "(t:%d) Pv: %f v:%d\n", t,val, v);
				//fprintf(stderr, "%f +%f +%f %f \n", log(1.0-d0), Pv[v], log(m0->A[v][u]+ZERO),  pdfL(m0, u, O[t]));
				if(val > maxPv){ maxPv = val; maxv = v; }
			}//v
			if(maxPj > maxPv){ 	// if switch
				Pu[u] = maxPj;
				nSu[u] = _copy_path(Sj[maxj], nSj[maxj], Su[u]);
				Su[u][nSu[u]] = t; nSu[u]++;
			}else{			// else
				Pu[u] = maxPv;
				nSu[u] = _copy_path(Sv[maxv], nSv[maxv], Su[u]);
			}
			//fprintf(stderr, "P(u:%d,t:%d) %f\n", u, t, Pu[u]);
			//_print_path(stderr, Su[u], nSu[u]);
		}//u


		///---  Pi(t)  ---///
		// regime-switch
		int maxv = _findMax(Pv, k0);
		for(i=0;i<k1;i++){
			// if switch
			double maxPv = Pv[maxv] + log(d0) + log(m1->pi[i]+ZERO)   + pdfL(m1, i, O[t]);
			//fprintf(stderr, "(t:%d) Pv: %f v:%d\n", t, maxPv, maxv);
			//fprintf(stderr, "%f +%f +%f +%f \n", Pv[maxv], + log(d0), + log(m1->pi[i]+ZERO),   + pdfL(m1, i, O[t]));
			// else
			double maxPj = -INF; int maxj;
			for(j=0;j<k1;j++){
				double val = log(1.0-d1) + Pj[j] + log(m1->A[j][i]+ZERO)  + pdfL(m1, i, O[t]);
				//fprintf(stderr, "(t:%d) Pj: %f j:%d\n", t,val, j);
				//fprintf(stderr, "%f +%f +%f %f \n", log(1.0-d1), Pj[j], + log(m1->A[j][i]+ZERO),  + pdfL(m1, i, O[t]));
				if(val > maxPj){ maxPj = val; maxj = j; }
			}//j
			if(maxPv > maxPj){ 	// if switch
				Pi[i] = maxPv;
				nSi[i] = _copy_path(Sv[maxv], nSv[maxv], Si[i]);
				Si[i][nSi[i]] = t; nSi[i]++;
			}else{			// else
				Pi[i] = maxPj;
				nSi[i] = _copy_path(Sj[maxj], nSj[maxj], Si[i]);
			}
			//fprintf(stderr, "P(i:%d,t:%d) %f\n", i, t, Pi[i]);
			//_print_path(stderr, Si[i], nSi[i]);
		}//i

		// swap Pu Pv; Pi Pj;
		dvec = Pu; Pu = Pv; Pv = dvec;
		dvec = Pi; Pi = Pj; Pj = dvec;
		// swap Su Sv; Si Sj;
		imat = Su; Su = Sv; Sv = imat;
		imat = Si; Si = Sj; Sj = imat;
		// swap nSu nSv; nSi nSj;
		ivec = nSu; nSu = nSv; nSv = ivec;
		ivec = nSi; nSi = nSj; nSj = ivec;
		//fprintf(stderr, "Pv:"); _printP(Pv, k0);fprintf(stderr, "Pj:");_printP(Pj, k1);

	}//t

	//fprintf(stderr, "Pv:"); _printP(Pv, k0);fprintf(stderr, "Pj:");_printP(Pj, k1);
	// check the paths
	int maxv = _findMax(Pv, k0);
	int maxj = _findMax(Pj, k1);
	//fprintf(stderr, "maxv:%d, maxj:%d\n", maxv, maxj);	
	int *path; int npath; int firstID;
	double lh = 0;
	if(Pv[maxv] > Pj[maxj]){
		path = Sv[maxv]; npath = nSv[maxv]; firstID=pow(-1,npath)*S0; lh = Pv[maxv];
	}else{
		path = Sj[maxj]; npath = nSj[maxj]; firstID=pow(-1,npath)*S1; lh = Pj[maxj];
	}
	//_print_path(stderr, path, npath); fprintf(stderr, "ID:%d\n");

	// add paths
	int curSt=st; int nxtSt;
	for(i=0;i<npath;i++){
		nxtSt = path[i];
		if(firstID*pow(-1,i)==S0)
			AddStEd(s0, curSt, nxtSt-curSt);
		else
			AddStEd(s1, curSt, nxtSt-curSt);
		curSt = nxtSt;	
	}//i
	if(firstID*pow(-1,npath)==S0)
		AddStEd(s0, curSt, st+len-curSt);
	else
		AddStEd(s1, curSt, st+len-curSt);

	double costC = -lh/log(2.0);

	return costC;
}



