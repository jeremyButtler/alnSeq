#include "strToNum.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){ /*main*/
  ulong maxValUL = defMaxUL();
  ulong tmpUL = 0;

  ulong uiCnt = 0;
  ulong numRoundsUL = 10000000;
  char uiToCnvtStr[23];
  char *tmpStr = 0;

  clock_t startClock;
  clock_t endClock;
  double timeUsedDbl;

  sprintf(uiToCnvtStr, "%lu", (ulong) 2947591200);
  strToULBase10(uiToCnvtStr, tmpUL);
  printf(
     "2 billon test (expect %s): %lu\n",
     uiToCnvtStr,
     tmpUL
  );

  sprintf(uiToCnvtStr, "%lu", (ulong) maxValUL);
  strToULBase10(uiToCnvtStr, tmpUL);
  printf("Max UL test (expect %lu): %lu\n",defMaxUL(),tmpUL);

  /*Overflow error; extra digit*/
  uiToCnvtStr[20] = '1';
  uiToCnvtStr[21] = '\0';
  tmpStr = strToULBase10(uiToCnvtStr, tmpUL);
  printf(
     "Overflow extra digit test (expect: %lu): %lu\n",
     defMaxUL(),
     tmpUL
  );
  printf("Extra digit should be 1: %s\n", tmpStr);

  /*Overflow error; 1 beyond max*/
  uiToCnvtStr[20] = '\0';
  ++uiToCnvtStr[19];
  tmpStr = strToULBase10(uiToCnvtStr, tmpUL);
  printf(
     "Overflow one dig off test (expect: %lu): %lu\n",
     defMaxUL() / 10,
     tmpUL
  );
  printf("Extra digit should be 6: %s\n", tmpStr);

  --uiToCnvtStr[19];

  startClock = clock();

  for(uiCnt = 0; uiCnt <= numRoundsUL; ++uiCnt)
  { /*Loop: Get number to convert*/
     sprintf(uiToCnvtStr, "%lu", uiCnt);
     strToULBase10(uiToCnvtStr, tmpUL);
     /*printf("%lu\n", tmpUL);*/
  } /*Loop: Get number to convert*/

  endClock = clock();
  timeUsedDbl =
     ((double) endClock - startClock) / CLOCKS_PER_SEC;

  printf("strToNum UL: %f\n", timeUsedDbl);

  startClock = clock();

  for(uiCnt = 0; uiCnt <= numRoundsUL; ++uiCnt)
  { /*Loop: Get number to convert*/
     sprintf(uiToCnvtStr, "%lu", uiCnt);
     tmpUL = strtoul(uiToCnvtStr, &tmpStr, 10);
     /*printf("%lu\n", tmpUL);*/
  } /*Loop: Get number to convert*/

  endClock = clock();
  timeUsedDbl =
     ((double) endClock - startClock) / CLOCKS_PER_SEC;

  printf("strtoul: %f\n", timeUsedDbl);

  exit(0);
} /*main*/
