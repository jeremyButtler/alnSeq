#include "alnMatrixStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 initAlnMatrixST:
'    - Sets all variables in matrixST to 0
'  o fun-02 freeAlnMatrixST
'    - Sets all variables in matrixST to 0
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o All variables in matrixST to be 0
\--------------------------------------------------------*/
void initAlnMatrixST(
  struct alnMatrixStruct *matrixST // Struct to initialize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-01: initAlnMatrixST
   '  - Sets all variables in matrixST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   matrixST->dirMatrixST = 0; 
   matrixST->refScoresST = 0;
   matrixST->lenRefScoresUL = 0;
   matrixST->queryScoresST = 0;
   matrixST->lenQueryScoresUL = 0;
   initScoresST(&matrixST->bestScoreST);

   return;
} // initAlnMatrixST

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o Frees alnMatrix and all of its held variables
\--------------------------------------------------------*/
void freeAlnMatrixST(
  struct alnMatrixStruct *matrixST // Struct to initialize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-0w TOC: Sec-01: freeAlnMatrixST
   '  - Sets all variables in matrixST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct scoresStruct *scoreST = 0;

   if(matrixST->dirMatrixST != 0)
     freeTwoBitAry(matrixST->dirMatrixST, 0, 0);
     // 0's to mark all elements on the heap

   if(matrixST->refScoresST != 0)
   { // if I need to free the reference scores
     scoreST = *matrixST->refScoresST;

     for(
       unsigned long ulScore = 0;
       ulScore < matrixST->lenRefScoresUL;
       ++ulScore
     ){ // for loop free all scores structs
       freeScoresST(scoreST, 0); // 0 for not on stack
       ++scoreST;
     } // for loop free all scores structs
   } // if I need to free the reference scores

   if(matrixST->queryScoresST != 0)
   { // if I need to free the query scores
     scoreST = *matrixST->queryScoresST;

     for(
       unsigned long ulScore = 0;
       ulScore < matrixST->lenQueryScoresUL;
       ++ulScore
     ){ // for loop free all scores structs
       freeScoresST(scoreST, 0); // 0 for not on stack
       ++scoreST;
     } // for loop free all scores structs
   } // if I need to free the query scores

   return;
} // freeAlnMatrixST
