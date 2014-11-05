#define YES 1
#define NO 0

#define FB 4*8 // float: 4*8;


double log_2(int x);
double log_s(int x);
double costHMM(int k, int d);

/// tool.c
void   setseed(void); 
double getrand(void);



int dFindMin(int *idx, double *value, int len);
int dFindMax(int *idx, double *value, int len);

int count_lines();
double getrusage_sec();
void notfin();
void error();
void SaveModel();

void LoadDB();
void NormSequence();
void ZnormSequence();

void PrintInput();
void PrintIdx();
void PrintSequence();
