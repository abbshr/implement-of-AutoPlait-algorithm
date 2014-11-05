/* 
 * plait.c --- written by Yasuko Matsubara 
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
#define DBG 0 //1
#define INFER_ITER_MIN 3  // #of min iter 
#define INFER_ITER_MAX 10 // #of max iter 
// minimum/maximum value for optimization
#define MINK  1 // minimum number of K
#define LM 0.1  // for sampling 
#define REGIME_R  0.03
#define SEGMENT_R 0.03 
//output option
#define PRINTHMM 1 
#define VITPATH 0 

PlaitWS *ws;
/// main function ///
void _plait();
void _regimeSplit();
void _cps();

/// find regimes & model estimation ///
void _regimeEst_aux();
void _estimateHMM_k();
void _estimateHMM();
double _findCentroid();
int _findOptSeedLen();
void _removeNoise();
void _selectLargest();

/// compute Vit & MDL ///
double _viterbi();
double _MDLtotal();
void _computeViterbiPath(Stack *s); //optional
void _computeLhMDL();
double _MDL();
/// output ///
void ReportPlait();
void _printSplit();
/// segment set --- get & release ///
SegBox *_getS();
void _releaseS();
/// memory allocation ///
void AllocPlait();



/// main function ///
void Plait(PlaitWS *pws){
	ws = pws;
	fprintf(stdout, "---------\n"); 
	fprintf(stdout, "r|m|Cost \n"); 
	fprintf(stdout, "---------\n"); 
	_plait();
}
/// main function ///

void _plait(){
	int k;
	double costT_s01;
	/// initialize Sx (X[0:m])
	SegBox *Sx = _getS("", "");
	Sx->costT= INF;
	AddStEd(Sx, 0, ws->x->m);		
	_estimateHMM(Sx);
	Push(Sx,&ws->C);
	while(1){
		ws->costT = _MDLtotal(&ws->Opt, &ws->C);
		Sx = Pop(&ws->C);
		if(Sx == NULL) break;
		/// s0 s1, create new segsets
		SegBox *s0 = _getS(Sx->label, "0");
		SegBox *s1 = _getS(Sx->label, "1");	
		/// try to split regime: Sx->(s0,s1)
		_regimeSplit(Sx, s0, s1);
		double costT_s01=s0->costT+s1->costT;
		/// split or not 
		if(costT_s01 + Sx->costT*REGIME_R < Sx->costT){
			Push(s0, &ws->C);
			Push(s1, &ws->C);
			_releaseS(Sx);
		}else{
			Push(Sx, &ws->Opt);
			_releaseS(s0); _releaseS(s1);
		}//if
	}//while
}
void _regimeSplit(SegBox *Sx, SegBox *s0, SegBox *s1){
	int seedlen= (int)ceil(ws->lmax*LM); 
	//int seedlen = _findOptSeedLen(Sx);
	// initialize HMM params
	_findCentroid(Sx, s0, s1, NSAMPLE, seedlen);
	_regimeEst_aux(Sx, s0, s1);
	if(s0->nSb==0 || s1->nSb==0)
		return;	
	/// final model estimation
	_estimateHMM(s0); _estimateHMM(s1);
	#if(DBG)
	_printSplit(stderr, Sx, s0, s1);
	#endif
}
void _regimeEst_aux(SegBox *Sx, SegBox *s0, SegBox *s1){
	int i;
	SegBox *opt0 = _getS("",""), *opt1 = _getS("","");
	for(i=0;i<INFER_ITER_MAX;i++){
		///--- Phase 1: Estimate parameters
		_selectLargest(s0); _selectLargest(s1);
		_estimateHMM(s0); _estimateHMM(s1);
		///--- Phase 2: Find cut-points	
		_cps(Sx, s0, s1, YES);
		if(s0->nSb==0 || s1->nSb==0) break; // avoid null inference
		/// if improving, update the opt seg-set
		double diff = (opt0->costT+opt1->costT) - (s0->costT+s1->costT);
		if(diff >0){ CopyStEd(s0, opt0); CopyStEd(s1, opt1); }
		// if not improving, then break iteration (efficient convergent)
		else if(i>=INFER_ITER_MIN)  break;  
	}//iter
	CopyStEd(opt0, s0); CopyStEd(opt1, s1);
	_releaseS(opt0); _releaseS(opt1);
	#if(DBG)
	fprintf(stderr, "=============== #iter:%d\n", i);
	#endif

}



void _selectLargest(SegBox *s){
	int id  = FindMaxStEd(s);
	int st  = s->Sb[id].st;
	int len = s->Sb[id].len;
	ResetStEd(s);
	AddStEd(s, st, len);
}

int _findMinDiff(SegBox *s0, SegBox *s1, double *diffp){
	double min = INF;
	int loc = -1;
	int i;
	for(i=0;i<s0->nSb;i++){
		int st  = s0->Sb[i].st;
		int len = s0->Sb[i].len;
		double costC0 = _viterbi(&s0->model, s0->delta, st, len, &ws->vit);
		double costC1 = _viterbi(&s1->model, s1->delta, st, len, &ws->vit);
		double diff = costC1-costC0;
		if(min > diff) {loc = i; min = diff;}
	}//i
	*diffp = min;	
	return loc;
}
double _scanMinDiff(SegBox *Sx, SegBox *s0, SegBox *s1){
	double diff;
	int loc0 = _findMinDiff(s0, s1, &diff);
	int loc1 = _findMinDiff(s1, s0, &diff);
	if(loc0==-1 || loc1==-1) {return INF;}
	SegBox *tmp0 = _getS("","");
	SegBox *tmp1 = _getS("","");
	AddStEd(tmp0, s0->Sb[loc0].st, s0->Sb[loc0].len);
	AddStEd(tmp1, s1->Sb[loc1].st, s1->Sb[loc1].len);
	_estimateHMM_k(tmp0, MINK);
	_estimateHMM_k(tmp1, MINK);
	double costC = CPSearch(Sx, tmp0, tmp1, ws);
	_releaseS(tmp0); _releaseS(tmp1);
	//fprintf(stderr, "%d:%d, %d:%d, %.0f\n", s0->Sb[loc0].st, s0->Sb[loc0].len, s1->Sb[loc1].st, s1->Sb[loc1].len, costC);
	return costC;
}

// remove too noisy small segments
void _removeNoise_aux(SegBox *Sx, SegBox *s0, SegBox *s1, double per){
	if(per==0) return;
	int i;
	int mprev=INF;
	double th = ws->costT*per;
	while(mprev > s0->nSb+s1->nSb){
		mprev = s0->nSb+s1->nSb;
		///--- find minimum segment ---///
		double diff0, diff1;
		int loc0 = _findMinDiff(s0, s1, &diff0);
		int loc1 = _findMinDiff(s1, s0, &diff1);
		double min; int id;
		if(diff0 < diff1){ min = diff0; id = 0; } else{ min = diff1; id = 1;}
		// check remove or not
		if(min<th){
			if(id==0){ AddStEd(s1, s0->Sb[loc0].st, s0->Sb[loc0].len); RemoveStEd(s0,loc0); } 
			else     { AddStEd(s0, s1->Sb[loc1].st, s1->Sb[loc1].len); RemoveStEd(s1,loc1); }
		}//if
	}//while
}

void _removeNoise(SegBox *Sx, SegBox *s0, SegBox *s1){ 
	if(s0->nSb<=1 && s1->nSb<=1) return;
	/// default pruning 
	double per = SEGMENT_R;
	_removeNoise_aux(Sx, s0, s1, per);
	double costC=_scanMinDiff(Sx, s0, s1);
	/// opt: optimal segset
	SegBox *opt0= _getS("",""); SegBox *opt1 = _getS("","");
	CopyStEd(s0, opt0); CopyStEd(s1, opt1);
	double prev=INF;
	/// find optimal pruning point
	while(per<=SEGMENT_R*10){
		if(costC>=INF) break;
		per*=2;
		_removeNoise_aux(Sx, s0, s1, per);
		if(s0->nSb<=1 || s1->nSb<=1) break;
		costC=_scanMinDiff(Sx, s0, s1);
		if(prev > costC){
			CopyStEd(s0, opt0); CopyStEd(s1, opt1);
		}else break;
		prev = costC;
	}//while
	//fprintf(stderr, "per:%f, cost: %.0f\n", per, costC);
	CopyStEd(opt0, s0); CopyStEd(opt1,s1);
	_releaseS(opt0);  _releaseS(opt1); 
}

void _cps(SegBox *Sx, SegBox *s0, SegBox *s1, int RM){
	CPSearch(Sx, s0, s1, ws);
	if(RM) _removeNoise(Sx, s0, s1);
	_computeLhMDL(s0); _computeLhMDL(s1);
}
double _viterbi(HMM *phmm, double delta, int st, int len, VITERBI *vit){
	double Lh = ViterbiL(phmm, len, ws->x->O+st, vit);
	if(delta<=0 || delta >= 1) error("not appropriate delta", "_computeLhMDL");
	Lh += log(delta);		// switch
	Lh += (len-1)*log(1.0-delta);	// else (stay)
	double costC = -Lh/log(2.0); 
	return costC;
}


void _estimateHMM_k(SegBox *s, int k){
	if(k<MINK) k=MINK; if(k>MAXK) k=MAXK;
	int wantKMEANS=1;
	int i,j;
	int n = s->nSb;
	/// if nSb is too big, use small size
	if(n > MAXBAUMN)
		n=MAXBAUMN;
	for(i=0;i<n;i++)
		Split(ws->x, s->Sb[i].st, s->Sb[i].len, &ws->x_tmp[i]);
	s->model.k = k;
	s->model.n = 0;
	ResetHMM(&s->model, k, ws->d);
	BaumWelch(&s->model, n, ws->x_tmp, &ws->baum, wantKMEANS);
	s->delta = (double)s->nSb/(double)s->len; //ws->lmax;
}

void _estimateHMM(SegBox *s){
	s->costT = INF;
	int k, optk;
	for(k=MINK; k<=MAXK; k++){
		double prev = s->costT;
		_estimateHMM_k(s, k);
		_computeLhMDL(s);
		if(s->costT > prev){ optk=k-1; break; }
	}//k
	if(optk<MINK) optk=MINK;
	if(optk>MAXK) optk=MAXK;
	_estimateHMM_k(s, optk);
	_computeLhMDL(s);
}

int _findOptSeedLen(SegBox *Sx){
	SegBox *s0= _getS("","");
	SegBox *s1= _getS("","");
	int iter1, iter2;
	double min=INF, prev=INF, pprev=INF;
	int len = FB;
	int minlen = Sx->Sb[FindMinStEd(Sx)].len;
	int maxlen = (int)ceil(ws->lmax*LM);
	int optlen= minlen;
	while((double)len <= maxlen && len <= minlen){
		int lm = Sx->len;
		double cost = _findCentroid(Sx, s0, s1, NSAMPLE, len);
		cost += len*log_s(MINK);
		#if(DBG)
		fprintf(stderr, "%d %f\n", len, cost);
		#endif
		//if(cost> prev) break; prev = cost; optlen = len;
		if(pprev < prev && pprev < cost) {optlen=len/(2*2); break;}
		pprev=prev; prev = cost;
		if(cost< min){min = cost; optlen = len;}//if
		/// next length:  current*2, i.e., 2^x; 
		len*=2;
	}//while
	#if(DBG)
	fprintf(stderr, "optimal minimum length: %d\n", optlen);
	#endif
	_releaseS(s0);
	_releaseS(s1);
	return optlen;	
}
double _findCentroid(SegBox *Sx, SegBox *s0, SegBox *s1, int nsamples, int seedlen){
	int iter1, iter2;
	double costMin = INF;
	/// keep best seeds
	int s0stB, s1stB, s0lenB, s1lenB; //best
	int s0stC, s1stC, s0lenC, s1lenC; //current
	/// make sample set
	UniformSet(Sx, seedlen, nsamples, &ws->U);
	/// start uniform sampling
	for(iter1=0;iter1<ws->U.nSb;iter1++)
	  for(iter2=iter1+1;iter2<ws->U.nSb;iter2++){
		UniformSampling(Sx, s0, s1, seedlen, iter1, iter2, &ws->U);
		if(s0->nSb==0 || s1->nSb==0) continue; // not sufficient
		// copy positions
		s0stC=s0->Sb[0].st; s0lenC=s0->Sb[0].len;
		s1stC=s1->Sb[0].st; s1lenC=s1->Sb[0].len;
		// estimate HMM
		_estimateHMM_k(s0, MINK);
		_estimateHMM_k(s1, MINK);
		// cut point search
		_cps(Sx, s0, s1, YES);
		if(s0->nSb==0 || s1->nSb==0) continue;
		if(costMin > s0->costT + s1->costT){
			// update best seeds
			costMin = s0->costT+s1->costT;
			s0stB=s0stC; s0lenB=s0lenC;
			s1stB=s1stC; s1lenB=s1lenC;
		}//if		
	}//iter1,2	
	if(costMin==INF){
		//fprintf(stderr, "couldn't find good seeds...\n");
		FixedSampling(Sx, s0, s1, seedlen);
		return INF;
	}
	ResetStEd(s0); ResetStEd(s1);
	AddStEd(s0, s0stB, s0lenB);
	AddStEd(s1, s1stB, s1lenB);
	return costMin;
}


/// compute Viterbi
void _copyVitPath(int *q, int st, int len, VITERBI *vit){
	int i;
	for(i=0;i<len;i++)
		q[i+st] = vit->q[i];
}
void _computeViterbiPath(Stack *C){
	int i,j;
	for(j=0;j<C->idx;j++){
		SegBox *s = C->s[j];
		for(i=0;i < s->nSb; i++){
			int st = s->Sb[i].st;
			int len = s->Sb[i].len;
			s->costC += _viterbi(&s->model, s->delta, st, len, &ws->vit);
			_copyVitPath(ws->q, st, len, &ws->vit);
		}//i
	}//j
	for(i=0;i<ws->lmax;i++)
		fprintf(stdout, "%d ", ws->q[i]);
  	fprintf(stdout, "\n");
}

/// compute MDL
void _computeLhMDL(SegBox *s){
	int i;
	if(s->nSb==0){
		s->costC = INF;
		s->costT = INF;
		return;
	}//if
	s->costC = 0;
	for(i=0;i < s->nSb; i++){
		int st = s->Sb[i].st;
		int len = s->Sb[i].len;
		s->costC += _viterbi(&s->model, s->delta, st, len, &ws->vit);
	}//i
	s->costT= _MDL(s);
}
double _MDL(SegBox *s){
	double costT = 0.0;
	int k = s->model.k;
	int d = s->model.d;
	int m = s->len;
	int i;
	double costC = s->costC;
	double costM = costHMM(k, d); 
	double costLen = 0.0;
	for(i=0;i<s->nSb;i++)
		costLen += log_2(s->Sb[i].len); 
	costLen += m * log_2(k);
	costT = costC + costLen + costM;
	return costT;
}
double _MDLtotal(Stack *Opt, Stack *C){
	int r = Opt->idx + C->idx;
	int m = mSegment(Opt) + mSegment(C);			
	double cost = MDLSegment(Opt) + MDLSegment(C);
	double costT = cost + log_s(r) + log_s(m) + m*log_2(r) + FB*r*r;
	fprintf(stdout, "%d %d %.0f \n", r, m, costT); 
	return costT;
}

/// output ///
void _printSplit(FILE *fp, SegBox *Sx, SegBox *s0, SegBox *s1){
	fprintf(fp, "[%s] %.0f(k:%d): {[%s] %.0f [%s] %.0f = %.0f(k:%d,%d)}\n",
	Sx->label, Sx->costT, Sx->model.k, 
	s0->label, s0->costT, 
	s1->label, s1->costT, 
	s0->costT+s1->costT, s0->model.k, s1->model.k);  
	/// segment
	fprintf(fp, "====\n");
	PrintStEd(fp, s0);
	fprintf(fp, "-------------\n");
	PrintStEd(fp, s1);
	fprintf(fp, "====\n");
}
void SavePlait(PlaitWS *pws, char *fdir){
    FILE *fp;
    int i;
    char filename[BUFSIZE];
    ws = pws;
    /// print segment positions
    for(i=0;i<ws->Opt.idx;i++){
    	sprintf(filename, "%ssegment.%d",fdir, i); 
    	if(( fp = fopen (filename, "w") ) == NULL )
    	  error("can not open:",filename);
	PrintStEd(fp, ws->Opt.s[i]);
    	fclose(fp);
    }//i
    /// print regime labels
    sprintf(filename, "%ssegment.labels",fdir); 
    if(( fp = fopen (filename, "w") ) == NULL )
      error("can not open:",filename);
    for(i=0;i<ws->Opt.idx;i++)
      fprintf(fp, "%d\t\t%s\t\t%.0f\t\t%d \n", i, ws->Opt.s[i]->label, ws->Opt.s[i]->costT, ws->Opt.s[i]->model.k);
    fclose(fp);
#if(PRINTHMM)
    /// print HMM params
    for(i=0;i<ws->Opt.idx;i++){
    	sprintf(filename, "%smodel.%d",fdir, i); 
    	if(( fp = fopen (filename, "w") ) == NULL )
    	  error("can not open:",filename);
	PrintHMM(fp, &ws->Opt.s[i]->model);
    	fclose(fp);
    }//i
#endif
#if(VITPATH)
   /// print vit path
   _computeViterbiPath(&ws->Opt);
#endif
}

/// segment set --- get & release ///
SegBox *_getS(char *parent, char *label){
	SegBox *s = Pop(&ws->S);
	if(s==NULL)
		error("too small maxc", "_getS()");
	ResetStEd(s);
	ResetHMM(&s->model, MAXK, ws->d);
	sprintf(s->label, "%s%s", parent, label);
	s->delta = 1.0/(double)ws->lmax; 
	return s;
}
void _releaseS(SegBox *s){
	Push(s, &ws->S);
}

/// memory alloc ///
void _allocSegBox(SegBox *s, int n){
  s->Sb = (SubS*)malloc(sizeof(SubS)*n);
  s->nSb = 0;
  s->optimal = 0;
  InitHMM(&s->model, MAXK, ws->d);
}
void AllocPlait(PlaitWS *pws){ 
  fprintf(stderr, "memory allocation...\n");
  ws = pws;
  ws->maxc   = (int) (log(ws->lmax) + 20); 
  ws->maxseg = (int) ws->lmax*0.1; 
  // alloc Viterbi
  AllocViterbi(&ws->vit, ws->maxseg, ws->lmax, MAXK);
  AllocViterbi(&ws->vit2, ws->maxseg, ws->lmax, MAXK);
  ws->q = (int*) ivector(0, ws->lmax);
  AllocCPS(&ws->cps, MAXK, ws->lmax);
  // alloc Baum
  int n = ws->maxseg;
  if(n > MAXBAUMN) n = MAXBAUMN;
  ws->x_tmp = (Input*)malloc(sizeof(Input)*n);
  AllocBaum(&ws->baum, n, ws->lmax, MAXK, ws->d);
  int i;
  // alloc segbox
  ws->s = (SegBox*)malloc(sizeof(SegBox)*ws->maxc);
  for(i=0;i<ws->maxc;i++)
	_allocSegBox(&ws->s[i], ws->maxseg);
  /// candidate stack C
  ws->C.s =(SegBox**)malloc(sizeof(SegBox*)*ws->maxc);
  ws->C.idx = 0;
  ws->Opt.s =(SegBox**)malloc(sizeof(SegBox*)*ws->maxc);
  ws->Opt.idx = 0;
  /// segment set storage
  ws->S.s =(SegBox**)malloc(sizeof(SegBox*)*ws->maxc);
  ws->S.idx = 0;
  for(i=0;i<ws->maxc;i++)
  	Push(&ws->s[i], &ws->S);
  /// for uniform sampling
  _allocSegBox(&ws->U, NSAMPLE);
}



