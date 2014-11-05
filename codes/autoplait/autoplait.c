/* 
 * autoplait.c --- written by Yasuko Matsubara 
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
#define DBG     1 
#define LMAX 10000000

int main(int argc,char *argv[]){
  double start,end,totaltime;
  if(argc != 4 + LSET)
	error("usage: cmd d fin, fout", argv[0]);
  setseed();
  PlaitWS ws;
  ws.lmax = LMAX; 

  //----get the "filename, dim, bin ..."----//
  ws.d = atoi(argv[1]);
  char* fin  = argv[2];
  char* fout = argv[3];
  if(LSET) ws.lmax = atoi(argv[4]);
  //------------------------------------//

  fprintf(stderr, "loading...\n");
  LoadDB(fin, 1, ws.lmax, ws.d, &ws.x);
  ws.lmax = ws.x->m;
  fprintf(stderr, "duration:  %d\n", ws.lmax);
  fprintf(stderr, "dimension: %d\n", ws.d);
  ZnormSequence(ws.x, 1, ws.d);
  AllocPlait(&ws);

  //---time check---//
  start = getrusage_sec();
  //---time check---//
  fprintf(stderr, "start autoplait...\n");
  Plait(&ws);
  //------------------//
  end = getrusage_sec();
  totaltime = end-start;
  //------------------//

  SavePlait(&ws, fout);
  printf("==================================\n");
  printf("duration: %d\n", ws.lmax);
  printf("search time: %.8f sec.\n",totaltime);
  printf("total patterns: %d \n", ws.Opt.idx);
  printf("total segments: %d \n", mSegment(&ws.Opt));
  printf("total cost: %.0f \n", ws.costT);
  printf("==================================\n");


}


