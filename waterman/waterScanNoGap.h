/*########################################################
# Name waterQryRefScan
# Use:
#  o Has the Waterman Smith query reference scan functions
#    that do not use gap extension penalties
# Libraries:
#   - "genWaterScanNoGap.h"            (No .c file)
#   o "../general/genScan.h"           (No .c file)
#   o "../general/genMath.h"           (No .c file)
#   o "../general/alnMatrixStruct.h"   (No .c file)
#   o "../general/twoBitArrays.h"      (No .c file)
#   o "../general/dataTypeShortHand.h" (No .c file)
#   - "../general/alnSetStruct.h"      (No .c file)
#   o "../general/base10StrToNum.h"    (No .c file)
#   o "../general/alnSeqDefaults.h"    (No .c file)
#   - "../general/seqStruct.h"         (No .c file)
# C Standard libraries:
#   o <stdlib.h>
#   o <stdio.h>
########################################################*/

#ifndef WATER_SCAN_NO_GAP_H
#define WATER_SCAN_NO_GAP_H

#include "genWaterScanNoGap.h"
#include "../general/alnSetStruct.h"
#include "../general/seqStruct.h"

/*-------------------------------------------------------\
| Fun-01 TOC: WaterScanNoGap
|  - Run a Waterman Smith alignment that saves the best
|    score for each query and reference base in the
|    alignment on the input sequences. This does not use
|    gap extension penalties
| Input:
|  - qryST:
|    o Pionter to seqStruct structure with query sequence 
|    o offsetUL; were to start query alignment, is index 0
|    o endAlnUL; were to stop query alignment, is index 0
|  - refST:
|    o Pionter to seqStruct with reference sequence 
|    o offsetUL; were to start alignment on ref, index 0
|    o endAlnUL; were to stop alignment on ref, is index 0
|  - alnSet:
|    o Point to alnSet structure with alignment settings
| Output:
|  - Returns:
|    o alnMatrixStruct with the direction matrix & scores
|    o 0 for memory allocation errors
\-------------------------------------------------------*/
static struct alnMatrix * WaterScanNoGap(
    struct seqStruct *qryST, /*query sequence and data*/
    struct seqStruct *refST, /*ref sequence and data*/
    struct alnSet *settings  /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: WaterScanNoGap
   '  - Run a Waterman Smith alignment on input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in initial negatives for ref
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   '  o fun-01 sec-05:
   '    - Set up for returing rap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01: Variable declerations
   ^  o fun-01 sec-01 sub-01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun-01 sec-01 sub-02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun-01 sec-01 sub-03:
   ^    - Directinol matrix variables
   ^  o fun-01 sec-01 sub-04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-01 Sub-01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   /*Get start & end of query and reference sequences*/
   char *refSeqStr = 0;
   char *qrySeqStr = 0;

   /*Find the length of the reference and query*/
   ulong lenQryUL = qryST->endAlnUL - qryST->offsetUL + 1;
   ulong lenRefUL = refST->endAlnUL - refST->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   ulong lenMatrixUL = (lenRefUL + 1) * (lenQryUL + 1);
     /*+1 for the gap column and row*/

   ulong ulRef = 0;
   ulong ulQry = 0;

   /*****************************************************\
   * Fun-01 Sec-01 Sub-02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   long delScoreL = 0;   /*Score for doing an deletion*/
   long nextSnpScoreL = 0;/*Score for the next match/snp*/

   /*Marks when to reset score buffer (every second row)*/
   long *scoreAryL = 0; /*matrix to use in alignment*/

   /*****************************************************\
   * Fun-01 Sec-01 Sub-03:
   *  - Directinol matrix variables
   \*****************************************************/

   struct alnMatrix *retMatrixST = 0;
   char *dirMatrix = 0;/*Direction matrix*/
   char *insDir = 0;   /*Direction above cell*/

   /*For recording alternative alignments*/
   ulong indexUL = 0;      /*Index in matrix*/
   ulong *indexAryUL=0;    /*Row of starting indexes*/
   ulong *oldIndexAryUL=0; /*Last round starting indexes*/
   ulong *swapPtrUL = 0;   /*For swapping ulongs*/

   /*For keeping track of best scores*/
   ulong lenAltRowUL = 0;

   long *refScoreAryL = 0;
   ulong *refIndexAryUL = 0;
   ulong *refEndIndexAryUL = 0;

   long *qryScoreAryL = 0;
   ulong *qryIndexAryUL = 0;
   ulong *qryEndIndexAryUL = 0;

   ulong ulScore = 0; /*For loop at end*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   ^  o fun-01 sec-02 sub-01:
   ^    - Allocate memory for the alignment
   ^  o fun-01 sec-02 sub-02:
   ^    - Allocate memory for alternative alignments
   ^  o fun-01 sec-02 sub-03:
   ^    - Get memory for keeping track of starting indexes
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Allocate memory for the alignment
   \*****************************************************/

   retMatrixST = malloc(sizeof(struct alnMatrix));
   if(retMatrixST == 0) return 0;
   initAlnMatrix(retMatrixST);

   dirMatrix = malloc((lenMatrixUL +1) * sizeof(char));

   if(dirMatrix == 0)
   { /*If: Memory error*/
      freeAlnMatrix(retMatrixST);
      return 0;
   } /*If: Memory error*/

   retMatrixST->dirMatrix = dirMatrix;
   scoreAryL = calloc((lenRefUL + 1), sizeof(long));
      /*+ 1 is for the indel column*/

   if(scoreAryL == 0)
   { /*If: I had a memory error*/
     freeAlnMatrix(retMatrixST);
     return 0;
   } /*If: I had a memory error*/

   retMatrixST->lenArraysUL = lenRefUL + lenQryUL;

   retMatrixST->lenRefUL = lenRefUL;
   retMatrixST->refOffsetUL = refST->offsetUL;

   retMatrixST->lenQryUL = lenQryUL;
   retMatrixST->qryOffsetUL = qryST->offsetUL;

   /*****************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Allocate memory for alternative alignments
   \*****************************************************/

   /*Length of each array (all same length)*/
   retMatrixST->lenArraysUL = lenRefUL + lenQryUL + 1;
   lenAltRowUL = retMatrixST->lenArraysUL;

   /*Set up the best scores array*/
   refScoreAryL =
      calloc(lenAltRowUL, sizeof(long));

   if(refScoreAryL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      return 0;
   } /*If had a memory error*/

   retMatrixST->scoreAryL = refScoreAryL;
   qryScoreAryL = refScoreAryL + lenRefUL;

   /*Set up the row of starting indexes*/
   refIndexAryUL =
      calloc(lenAltRowUL, sizeof(ulong));

   if(refIndexAryUL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      return 0;
   } /*If had a memory error*/

   retMatrixST->startIndexAryUL = refIndexAryUL;
   qryIndexAryUL = refIndexAryUL + lenRefUL;

   /*Set up the row of ending indexes*/
   refEndIndexAryUL =
      calloc((lenAltRowUL), sizeof(ulong));

   if(refEndIndexAryUL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      return 0;
   } /*If had a memory error*/

   retMatrixST->endIndexAryUL = refEndIndexAryUL;
   qryEndIndexAryUL = refEndIndexAryUL + lenRefUL;

   /*****************************************************\
   * Fun-01 Sec-02 Sub-03:
   *  - Get memory for keeping track of starting indexes
   \*****************************************************/

   /*Set up the first row of starting indexes*/
   indexAryUL = malloc((lenRefUL + 1) * sizeof(ulong));

   if(indexAryUL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      return 0;
   } /*If had a memory error*/

   /*Set up the second row of indexs (so have two rows)*/
   oldIndexAryUL = malloc((lenRefUL + 1) * sizeof(ulong));

   if(oldIndexAryUL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      free(indexAryUL);
      return 0;
   } /*If had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   insDir = dirMatrix;

   for(indexUL = 0; indexUL <= lenRefUL; ++indexUL)
   { /*Loop: Set up the gap row*/
      dirMatrix[indexUL] = defMvStop;
      indexAryUL[indexUL] = indexUL;
   } /*Loop: Set up the gap row*/

   /*Scores handled by calloc*/
   dirMatrix += indexUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Final set up for socring loops
   ^  o fun-01 sec-04 sub-02:
   ^    - Fill out the matrix
   ^  o fun-01 sec-04 sub-03:
   ^    - Find scores for second half of reference
   ^  o fun-01 sec-04 sub-04:
   ^    - Find the next matches score
   ^  o fun-01 sec-04 sub-05:
   ^    - Find the best score for last base of each row
   ^  o fun-01 sec-04 sub-06:
   ^    - Prepare for the next round (row)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Final set up for scoring loop
   \*****************************************************/

   /*Setup the index arrays for the first scoring round*/
   swapPtrUL = indexAryUL;
   indexAryUL = oldIndexAryUL;
   oldIndexAryUL = swapPtrUL;

   /*Set up scores for the first scoring round*/
   nextSnpScoreL = scoreAryL[0];
   delScoreL = 0;
   indexAryUL[0] = indexUL;
   dirMatrix[0] = defMvStop;

   /*Position sequences for scoring*/
   ++indexUL;
   refSeqStr = refST->seqCStr + refST->offsetUL - 1;
   qrySeqStr = qryST->seqCStr + qryST->offsetUL;

   /*****************************************************\
   * Fun-01 Sec-04 Sub-02:
   *  - Fill out the matrix
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      ulQry = 0;
      ulQry < lenQryUL;
      ++ulQry
   ){ /*Loop: compare query base against all ref bases*/

      for(
         ulRef = 1;
         ulRef < (lenRefUL) / 2;
         ++ulRef
      ){ /* Loop: Check the 1st half of reference bases*/

         waterMatrixScanMaxScoreNoGap(
            refSeqStr[ulRef],
            qrySeqStr[ulQry],
            scoreAryL[ulRef],
            dirMatrix[ulRef],
            insDir[ulRef],
            indexAryUL[ulRef],
            nextSnpScoreL,
            delScoreL,
            oldIndexAryUL[ulRef - 1],
            oldIndexAryUL[ulRef],
            indexAryUL[ulRef - 1],
            indexUL,
            settings
         );

         scanIfKeepScoreRef(
            scoreAryL[ulRef],
            dirMatrix[ulRef],
            indexAryUL[ulRef],
            indexUL,
            refScoreAryL[ulRef],
            refIndexAryUL[ulRef],
            refEndIndexAryUL[ulRef],
            qryScoreAryL[ulQry],
            qryIndexAryUL[ulQry],
            qryEndIndexAryUL[ulQry]
         );

         ++indexUL;
      } /* Loop: Check the 1st half of reference bases*/

      /**************************************************\
      * Fun-01 Sec-04 Sub-03:
      *  - Find scores for second half of the reference
      \**************************************************/

      for(
         ulRef = ulRef;
         ulRef < lenRefUL;
         ++ulRef
      ){ /* Loop: Check last half of reference bases*/
         waterMatrixScanMaxScoreNoGap(
            refSeqStr[ulRef],
            qrySeqStr[ulQry],
            scoreAryL[ulRef],
            dirMatrix[ulRef],
            insDir[ulRef],
            indexAryUL[ulRef],
            nextSnpScoreL,
            delScoreL,
            oldIndexAryUL[ulRef - 1],
            oldIndexAryUL[ulRef],
            indexAryUL[ulRef - 1],
            indexUL,
            settings
         );

         scanIfKeepScoreQry(
            scoreAryL[ulRef],
            dirMatrix[ulRef],
            indexAryUL[ulRef],
            indexUL,
            refScoreAryL[ulRef],
            refIndexAryUL[ulRef],
            refEndIndexAryUL[ulRef],
            qryScoreAryL[ulQry],
            qryIndexAryUL[ulQry],
            qryEndIndexAryUL[ulQry]
         );

         ++indexUL;
      } /* Loop: Check last half of reference bases*/

      /***************************************************\
      * Fun-01 Sec-04 Sub-04:
      *  - Find the best score for last base of each row
      \***************************************************/

      waterMatrixScanMaxEndRowNoGap(
         refSeqStr[ulRef],
         qrySeqStr[ulQry],
         scoreAryL[ulRef],
         dirMatrix[ulRef],
         insDir[ulRef],
         indexAryUL[ulRef],
         nextSnpScoreL,
         delScoreL,
         oldIndexAryUL[ulRef - 1],
         oldIndexAryUL[ulRef],
         indexAryUL[ulRef - 1],
         indexUL,
         settings
      );

     scanIfKeepScoreQry(
        scoreAryL[ulRef],
        dirMatrix[ulRef],
        indexAryUL[ulRef],
        indexUL,
        refScoreAryL[ulRef],
        refIndexAryUL[ulRef],
        refEndIndexAryUL[ulRef],
        qryScoreAryL[ulQry],
        qryIndexAryUL[ulQry],
        qryEndIndexAryUL[ulQry]
     );

     ++indexUL;

     /***************************************************\
     *  Fun-01 Sec-04 Sub-05:
     *   - Prepare for the next round
     \***************************************************/

      /*Set up for next row of the matirx*/
      insDir = dirMatrix;
      dirMatrix += ulRef + 1;

      /*Get scores set up*/
      nextSnpScoreL = scoreAryL[0];
      dirMatrix[0] = defMvStop;
      scoreAryL[0] += settings->gapOpenC;
      scoreAryL[0] &= -(scoreAryL[0] > 0);

      delScoreL = scoreAryL[0] + settings->gapOpenC;
      delScoreL &= -(delScoreL > 0);

      /*Swap index arrays so the current is last*/
      swapPtrUL = indexAryUL;
      indexAryUL = oldIndexAryUL;
      oldIndexAryUL = swapPtrUL;

      ++indexUL; /*Set index for the next base pair*/
     /*At this piont insDir is on the second base*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   ^  o fun-01 sec-05 sub-01:
   ^    - clean up
   ^  o fun-01 sec-05 sub-02:
   ^    - find the best score
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-05 Sub-01:
   *  - clean up
   \*****************************************************/

   /*Move back to the lower right conor cell
   ` This is not needed, but is nice.
   */
   insDir[ulRef + 1] = defMvStop;

   free(scoreAryL);
   free(indexAryUL);
   free(oldIndexAryUL);

   scoreAryL = 0;
   indexAryUL = 0;
   oldIndexAryUL = 0;

   /*****************************************************\
   * Fun-01 Sec-05 Sub-01:
   *  - Find the best score
   \*****************************************************/

   scoreAryL = retMatrixST->scoreAryL;
   indexAryUL = retMatrixST->startIndexAryUL;
   oldIndexAryUL = retMatrixST->endIndexAryUL;

   for(
      ulScore = 0;
      ulScore < retMatrixST->lenArraysUL;
     ++ulScore
   ){ /*Loop: Find the highest score*/
      if(scoreAryL[ulScore] > retMatrixST->bestScoreL)
      { /*If I have a new best score*/
         retMatrixST->bestScoreL = scoreAryL[ulScore];

         retMatrixST->bestStartIndexUL =
            indexAryUL[ulScore];

         retMatrixST->bestEndIndexUL =
            oldIndexAryUL[ulScore];
      } /*If I have a new best score*/
   } /*Loop: Find the highest score*/

   return retMatrixST;
} /*WaterScanNoGap*/

#endif
