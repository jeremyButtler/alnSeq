#include "strToNum.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){ /*main*/
  uint maxValUI = defMaxUI();
  uint tmpUI = 0;
  ulong tmpUL = 0;

  uint uiCnt = 0;
  uint numRoundsUI = 10000000;
  char uiToCnvtStr[20];
  char *tmpStr = 0;

  clock_t startClock;
  clock_t endClock;
  double timeUsedDbl;

  sprintf(uiToCnvtStr, "%lu", (ulong) 2947591200);
  strToUIBase10(uiToCnvtStr, tmpUI);
  printf(
     "2 billon test (expect %s): %u\n",
     uiToCnvtStr,
     tmpUI
  );

  sprintf(uiToCnvtStr, "%lu", (ulong) maxValUI);
  strToUIBase10(uiToCnvtStr, tmpUI);
  printf("Max UI test (expect %u): %u\n",defMaxUI(),tmpUI);

  /*Overflow error; extra digit*/
  uiToCnvtStr[10] = '1';
  uiToCnvtStr[11] = '\0';
  tmpStr = strToUIBase10(uiToCnvtStr, tmpUI);
  printf(
     "Overflow extra digit test (expect: %u): %u\n",
     defMaxUI(),
     tmpUI
  );
  printf("Extra digit should be 1: %s\n", tmpStr);

  /*Overflow error; 1 beyond max*/
  uiToCnvtStr[10] = '\0';
  ++uiToCnvtStr[9];
  tmpStr = strToUIBase10(uiToCnvtStr, tmpUI);
  printf(
     "Overflow one dig off test (expect: %u): %u\n",
     defMaxUI() / 10,
     tmpUI
  );
  printf("Extra digit should be 6: %s\n", tmpStr);

  --uiToCnvtStr[9];

  startClock = clock();

  for(uiCnt = 0; uiCnt <= numRoundsUI; ++uiCnt)
  { /*Loop: Get number to convert*/
     sprintf(uiToCnvtStr, "%u", uiCnt);
     strToUIBase10(uiToCnvtStr, tmpUI);
     /*printf("%u\n", tmpUI);*/
  } /*Loop: Get number to convert*/

  endClock = clock();
  timeUsedDbl =
     ((double) endClock - startClock) / CLOCKS_PER_SEC;

  printf("strToNum UI: %f\n", timeUsedDbl);

  startClock = clock();

  for(uiCnt = 0; uiCnt <= numRoundsUI; ++uiCnt)
  { /*Loop: Get number to convert*/
     sprintf(uiToCnvtStr, "%u", uiCnt);
     tmpUL = strtoul(uiToCnvtStr, &tmpStr, 10);
     /*printf("%u\n", tmpUI);*/
  } /*Loop: Get number to convert*/

  endClock = clock();
  timeUsedDbl =
     ((double) endClock - startClock) / CLOCKS_PER_SEC;

  printf("strtoul: %f\n", timeUsedDbl);

  exit(0);
} /*main*/
