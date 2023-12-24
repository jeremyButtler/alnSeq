/*########################################################
# Name needleTwoBitNoGap
# Use:
#  o Holds functions for doing a pairwise Needleman Wunsch
#    alignment. This uses two bit arrays, but not gap
#    extension penalites
# Libraries:
#   - "genNeedleNoGap.h"               (No .c file)
#   o "../general/genAln.h"            (No .c file)
#   o "../general/genMath.h"           (No .c file)
#   - "../general/alnMatrixStruct.h"   (No .c file)
#   o "../general/twoBitArrays.h"      (No .c file)
#   o "../general/dataTypeShortHand.h" (No .c file)
#   - "../general/alnSetStruct.h"      (No .c file)
#   o "../general/base10StrToNum.h"    (No .c file)
#   o "../general/alnSeqDefaults.h"    (No .c file)
#   - "../general/seqStruct.h"         (No .c file)
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  - fun-01 NeedleManWunschAln:
'    o Perform a Needleman-Wunsch alignment on the two
'      input sequences with two bit arrays and no gap
'      extension penalties
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef NEEDLE_TWO_BIT_NO_GAP_H
#define NEEDLE_TWO_BIT_NO_GAP_H

#include "genNeedleNoGap.h"
#include "../general/alnMatrixStruct.h"
#include "../general/seqStruct.h"
#include "../general/alnSetStruct.h"

/*-------------------------------------------------------\
| Fun-01: NeedleTwoBitNoGap
|  - Perform a Needleman-Wunsch alignment on the two input
|    sequences with two bit arrays and without gap
|    extension penalties
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
static struct alnMatrixTwoBit * NeedleTwoBitNoGap(
    struct seqStruct *qryST, /*query sequence and data*/
    struct seqStruct *refST, /*ref sequence and data*/
    struct alnSet *settings  /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: NeedleTwoBitNoGap
   '  - Perform a Needleman-Wunsch alignment on the two
   '    input sequences using two bit arrays and no gap
   '    extension penalties
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in the initial negatives for the reference
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Get start & end of query & reference sequences*/
   char *refSeqStr = refST->seqCStr + refST->offsetUL;
   char *qrySeqStr = qryST->seqCStr + qryST->offsetUL;

   /*Find the length of the reference and query. The +1
   ` is to account for offsetUL being index 0
   */
   ulong lenQryUL = qryST->endAlnUL - qryST->offsetUL + 1;
   ulong lenRefUL = refST->endAlnUL - refST->offsetUL + 1;
   ulong lenMatrixUL = (lenRefUL + 1) * (lenQryUL + 1);

   /*Variables for loops*/
   ulong ulRef = 0;
   ulong ulQry = 0;

   /*Scoring variables*/
   long nextSnpScoreL = 0;   /*Score for doing an match/snp*/
   long delScoreL = 0;   /*Score for doing an deletion*/

   long *scoreAryL = 0; /*matrix to use in alignment*/

   /*Direction matrix (1 cell holds a single direction)*/
   struct alnMatrixTwoBit *retMatrixST = 0;
   struct twoBitAry *dirMatrix = 0;  /*Direction matrix*/
   struct twoBitAry insDir;         /*Direction above*/

   uchar dirUC = 0;    /*holds choosen direction*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retMatrixST = malloc(sizeof(struct alnMatrixTwoBit));
   if(retMatrixST == 0) return 0;
   initAlnMatrixTwoBit(retMatrixST);

   dirMatrix = makeTwoBit((lenMatrixUL+1), 0);

   if(dirMatrix == 0)
   { /*If: Memory error*/
      freeAlnMatrixTwoBit(retMatrixST);
      return 0;
   } /*If: Memory error*/

   retMatrixST->dirMatrix = dirMatrix;

   scoreAryL = malloc((lenRefUL + 1) * sizeof(long));
      /*+ 1 is for the indel column*/

   if(scoreAryL == 0)
   { /*If: I had a memory error*/
     freeAlnMatrix(retMatrixST);
     return 0;
   } /*If: I had a memory error*/

   retMatrixST->lenRefUL = lenRefUL;
   retMatrixST->refOffsetUL = refST->offsetUL;

   retMatrixST->lenQryUL = lenQryUL;
   retMatrixST->qryOffsetUL = qryST->offsetUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in the initial negatives for the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Find the first two insertions (non-standard)*/
   cpTwoBitPos(dirMatrix, &insDir);
   changeTwoBitElm(dirMatrix, defMvStop);
   twoBitMvToNextElm(dirMatrix);
   changeTwoBitElm(dirMatrix, defMvDel);
   twoBitMvToNextElm(dirMatrix);

   scoreAryL[0] = 0;
   scoreAryL[1] = settings->gapOpenC;

   for(ulRef = 2; ulRef <= lenRefUL; ++ulRef)
   { /*loop; till have initalized the first row*/
     changeTwoBitElm(dirMatrix, defMvDel);
     twoBitMvToNextElm(dirMatrix);

     scoreAryL[ulRef] = 
        scoreAryL[ulRef - 1] + settings->gapOpenC;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Final preperation before scoring
   ^  o fun-01 sec-04 sub-02:
   ^    - Get scores for each row
   ^  o fun-01 sec-04 sub-03:
   ^    - Get the last score at the end of each row
   ^  o fun-01 sec-04 sub-04:
   ^    - Prepare for scoring the next row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Final preperation before scoreing
   \*****************************************************/

   refSeqStr = refST->seqCStr + refST->offsetUL - 1;
   qrySeqStr = qryST->seqCStr + qryST->offsetUL;

   nextSnpScoreL = 0;

   /*Fill in the current indel column for this row*/
   changeTwoBitElm(dirMatrix, defMvIns);
   twoBitMvToNextElm(dirMatrix);
   twoBitMvToNextElm(&insDir);

   scoreAryL[0] = settings->gapOpenC;
   delScoreL = settings->gapOpenC;

   /*****************************************************\
   * Fun-01 Sec-04 Sub-02:
   *  - Get scores for each row
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      ulQry = 0;
      ulQry < qryST->endAlnUL - qryST->offsetUL + 1;
      ++ulQry
   ){ /*loop; fill the direction matrix with socres*/

       for(
          ulRef = 1;
          ulRef <= refST->endAlnUL - refST->offsetUL;
          ++ulRef
       ){ /*Loop: compare 1 query to all reference bases*/
          needleMaxScoreNoGap(
             refSeqStr[ulRef],
             qrySeqStr[ulQry],
             scoreAryL[ulRef],
             dirUC,
             getTwoBitElm(&insDir),
             nextSnpScoreL,
             delScoreL,
             settings 
          );

          changeTwoBitElm(dirMatrix, dirUC);
          twoBitMvToNextElm(dirMatrix);
          twoBitMvToNextElm(&insDir);
       } /*Loop: compare 1 query to all reference bases*/

       /*************************************************\
       * Fun-01 Sec-04 Sub-03:
       *  - Get scores for the end of each row
       **************************************************/

       needleMaxEndRowScoreNoGap(
          refSeqStr[ulRef],
          qrySeqStr[ulQry],
          scoreAryL[ulRef],
          dirUC,
          getTwoBitElm(&insDir),
          nextSnpScoreL,
          delScoreL,
          settings 
       );

       changeTwoBitElm(dirMatrix, dirUC);
       twoBitMvToNextElm(dirMatrix);
       twoBitMvToNextElm(&insDir);

       /*************************************************\
       * Fun-01 Sec-04 Sub-4:
       *  - Set up for scoring the next row
       \*************************************************/

       changeTwoBitElm(dirMatrix, defMvIns);
       twoBitMvToNextElm(dirMatrix);
       twoBitMvToNextElm(&insDir);

       nextSnpScoreL = scoreAryL[0];
       scoreAryL[0] += settings->gapOpenC;
       delScoreL = scoreAryL[0] + settings->gapOpenC;
   } /*loop; fill the direction matrix with socres*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Set the best score to the cornor right cell*/
   retMatrixST->bestScoreL = scoreAryL[lenRefUL];
   
   twoBitMvBackOneElm(dirMatrix);
   changeTwoBitElm(dirMatrix, defMvStop);
   twoBitMvBackOneElm(dirMatrix);

   retMatrixST->bestEndIndexUL =
      twoBitGetIndex(dirMatrix);

   /*Clean UP*/
   free(scoreAryL);
   scoreAryL = 0;

   return retMatrixST;
} /*NeeldeTwoBit*/

#endif
