/* 
 * segbox.c --- written by Yasuko Matsubara 
 * SegBox and Stack
 *
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

void PrintStEd(FILE *fp, SegBox *s){
	int i;
	for(i=0;i<s->nSb;i++)
		fprintf(fp, "%d %d \n", 
			s->Sb[i].st, 
			s->Sb[i].st+s->Sb[i].len-1);
}
void ResetStEd(SegBox *s){ 
	s->nSb = 0;
	s->len = 0;
	s->costC = INF;
	s->costT = INF;
}
void RemoveStEd(SegBox *s, int id){
	if(id >= s->nSb) { fprintf(stderr, "%d: ", id); error("too large id","RemoveStEd"); }
	int i;
	int len = s->Sb[id].len;
	for(i=id;i<s->nSb-1;i++){
		s->Sb[i].st=s->Sb[i+1].st;
		s->Sb[i].len=s->Sb[i+1].len;
	}//i
	s->len-=len;
	s->nSb--;
}

void _rm_overlap(SegBox *s){
	int i;
	int curr=INF;
	// check overlapped segments
	while(curr > s->nSb){
		curr = s->nSb;
		for(i=0;i<s->nSb-1;i++){
			int st0 = s->Sb[i].st;
			int ed0 = st0 + s->Sb[i].len;
			int st1 = s->Sb[i+1].st;
			int ed1 = st1 + s->Sb[i+1].len;
			int ed;
			if(ed0 > ed1) ed = ed0; else ed = ed1;
			// if overlapped
			if(ed0+1 >= st1){
				RemoveStEd(s, i+1);
				s->Sb[i].len = ed-st0;
			}//if
		}//i
	}//while
}

void AddStEd(SegBox *s, int st, int len){
	if(len<=0) return; //error("len<=0", "addStEd"); //return;
	if(st<0)   st=0; //error("st < 0", "addStEd"); //st=0;
	// find loc
	int i; int loc = -1;
	for(i=0;i<s->nSb;i++){
		if(s->Sb[i].st > st)
			break;
	}//i
	loc = i;
	// add new segment
	for(i=s->nSb-1;i>=loc;i--){
		s->Sb[i+1].st  = s->Sb[i].st;
		s->Sb[i+1].len = s->Sb[i].len;
	}//i
	s->Sb[loc].st = st;
	s->Sb[loc].len = len;
	s->nSb++;
	// remove overlap
	_rm_overlap(s);
	// length sum
	s->len = 0;
	for(i=0;i<s->nSb;i++)
		s->len += s->Sb[i].len;
}

/// allow overlapped segments
void AddStEd_ex(SegBox *s, int st, int len){
	s->len += len;
	int next = s->nSb;
	s->Sb[next].st = st; 
	s->Sb[next].len = len; 
	s->nSb++;
}

int Push(SegBox *s, Stack *C){
	C->s[C->idx] = s;
	C->idx++;
	return C->idx;
}
SegBox *Pop(Stack *C){
        if(C->idx==0)
                return NULL;
        C->idx--;
        return C->s[C->idx];
}
double MDLSegment(Stack *C){
	int i;
	double c = 0;
	for(i=0;i<C->idx;i++)
		c += C->s[i]->costT;
	return c;                    
}
int mSegment(Stack *C){
	int i;
	int m = 0;
	for(i=0;i<C->idx;i++)
		m += C->s[i]->nSb;
	return m;                    
}
void CopyStEd(SegBox *from, SegBox *to){
	int i;
	for(i=0;i<from->nSb;i++){
		to->Sb[i].st = from->Sb[i].st;
		to->Sb[i].len = from->Sb[i].len;
	}
	to->nSb   = from->nSb;
	to->len   = from->len;
	to->costT = from->costT;
	to->costC = from->costC;
	to->delta = from->delta;
	strcpy(to->label, from->label);
}

int IsSameStEd(SegBox *a, SegBox *b){
	int i;
	if(a->nSb != b->nSb)
		return -1;
	for(i=0;i<a->nSb;i++){
		if(a->Sb[i].st !=b->Sb[i].st)
			return -1;
		if(a->Sb[i].len !=b->Sb[i].len)
			return -1;
	}//i
	return 1;
}

int FindMinStEd(SegBox *s){
        int i,j;
        int id=-1; int min=INF;
        for(i=0;i<s->nSb;i++)
		if(s->Sb[i].len < min){ 
			id=i; min=s->Sb[i].len; 
		}
	return id;
}
int FindMaxStEd(SegBox *s){
        int i,j;
        int id=-1; int max=-INF;
        for(i=0;i<s->nSb;i++)
		if(s->Sb[i].len > max){ 
			id=i; max=s->Sb[i].len; 
		}
	return id;
}

void FixedSampling(SegBox *Sx, SegBox *s0, SegBox *s1, int len){
	int loc, r;
	// init segments
	ResetStEd(s0);
	ResetStEd(s1);
	// segment s0
	loc = 0%Sx->nSb;
	r = Sx->Sb[loc].st; 
	AddStEd(s0, r, len);
	// segment s1
	loc = 1%Sx->nSb; 
	r = Sx->Sb[loc].st+(int)(Sx->Sb[loc].len/2); 
	AddStEd(s1, r, len);
}
void RandomSampling(SegBox *Sx, SegBox *s0, SegBox *s1, int len){
	int loc, r;
	// init segments
	ResetStEd(s0);
	ResetStEd(s1);
	// segment s0
	loc = (int) floor(getrand()*Sx->nSb);
	r = (int) floor(Sx->Sb[loc].st + getrand()*(Sx->Sb[loc].len-len));
	AddStEd(s0, r, len);
	// segment s1
	loc = (int) floor(getrand()*Sx->nSb);
	r = (int) floor(Sx->Sb[loc].st + getrand()*(Sx->Sb[loc].len-len));
	AddStEd(s1, r, len);
}
void UniformSet(SegBox *Sx, int len, int trial, SegBox *U){
	int i,j;
	int slideW = (int)ceil((Sx->len-len)/trial);
	// create uniform blocks
	ResetStEd(U);	
	for(i=0;i<Sx->nSb;i++){
		if(U->nSb >= trial) return;
		int st=Sx->Sb[i].st;
		int ed=st+Sx->Sb[i].len;
		for(j=0;j<trial;j++){
			int next = st+j*slideW;
			if(next+len>ed){
				int st=ed-len; if(st<0)st=0;
				AddStEd_ex(U, st, len);
				break;	
			}	
			AddStEd_ex(U, next, len);
		}//j
	}//i
}
void UniformSampling(SegBox *Sx, SegBox *s0, SegBox *s1, int len, int n1, int n2, SegBox *U){
	int i,j;
	// init segments
	ResetStEd(s0);
	ResetStEd(s1);
	i = (int) (n1 % U->nSb);
	j = (int) (n2 % U->nSb);
	int st0 = U->Sb[i].st;
	int st1 = U->Sb[j].st;
	/// if overlapped, then, ignore
	if(abs(st0-st1)<len) return;
	AddStEd(s0, st0, len);
	AddStEd(s1, st1, len);
}

