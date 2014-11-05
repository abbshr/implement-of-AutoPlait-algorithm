/* 
 * autoplait.h --- written by Yasuko Matsubara 
 */

#define DEFAULT 1 
#if(DEFAULT)
 #define MAXK     16  
 #define LSET     NO  
 #define MAXBAUMN 3
 #define NSAMPLE  10 // #of sampling 
#else
 #define MAXK     4
 #define LSET NO  //YES
 #define MAXBAUMN 1
 #define NSAMPLE  5
#endif

typedef struct {
  int st;
  int len;
}SubS; //SubSequence

typedef struct {
	SubS *Sb;
	int nSb;
	int len;
	double costT;   // total cost
	double costC;   // coding cost
	int optimal;
	char label[BUFSIZE];
	HMM model;
	double delta;
}SegBox; //Segment

typedef struct {
  SegBox **s;
  int idx;
}Stack; //SegBox List

typedef struct {
	double *Pu, *Pv, *Pi, *Pj;
	int **Su, **Sv, **Si, **Sj;
	int *nSu, *nSv, *nSi, *nSj;
}CPS;  //CutPointSearch

typedef struct {
  Input *x;
  int maxc;
  int maxseg;
  int d;
  int lmax;
  // cut point search
  CPS cps;
  // baum & viterbi
  BAUM baum;
  int *q; //q[m] ... Viterbi path
  VITERBI vit;
  VITERBI vit2;
  Input *x_tmp;
  SegBox *s;
  // candidate stack
  Stack C;
  Stack Opt;
  Stack S;
  double costT;
  SegBox U; // uniform sampling
}PlaitWS; // autoplait workspace

/// plait.c
void PrintParams();
void AllocPlait();
void Plait();
void SavePlait();

/// segbox.c
void PrintStEd(FILE *fp, SegBox *s);
void CopyStEd(SegBox *from, SegBox *to);
void ResetStEd(SegBox *s);
void AddStEd(SegBox *s, int st, int len);
void AddStEd_ex(SegBox *s, int st, int len);
void RemoveStEd(SegBox *s, int id);
int Push(SegBox *s, Stack *C);
SegBox *Pop(Stack *C);
int mSegment(Stack *C);
double MDLSegment(Stack *C);
int IsSameStEd(SegBox *a, SegBox *b);
int FindMaxStEd();
int FindMinStEd();
void FixedSampling();
void RandomSampling();
void UniformSet();
void UniformSampling();

/// cps.c
double CPSearch();
void AllocCPS();
void FreeCPS();

