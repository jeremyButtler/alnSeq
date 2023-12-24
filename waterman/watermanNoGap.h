/*########################################################
# Name watermanNoGap
# Use:
#  o Holds functions doing a Waterman-Smith pairwise
#    alignments.
# Includes:
#   - "genWaterNoGap.h"                (No .c file)
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
#   o <stdio.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - Has includes and definede variables
'  o fun-01 WatermanSmithAln:
'    - Perform a Waterman Smith alignment on input
'      sequences
'  o fun-02 printAltWaterAlns:
'    - Prints out the best aligment and the saved
'       alterantive alignments  (best alignment for each
'       base) to a file
'  o fun-03 updateDirScoreWaterSingle:
'    - Picks the best score and direction for the current
'      base pairs being compared in a Waterman-Smith
'      alignment
'  o fun-04 keepAltScore:
'    - Detemines if keeping a score for a row or column
'      in a query ref scan
'  o fun-05 addBestBaseScore:
'    - Adds a score and index to the kept scores list
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Has includes and definede variables
\-------------------------------------------------------*/

#ifndef WATERMAN_NO_GAP_H
#define WATERMAN_NO_GAP_H

#include "genWaterNoGap.h"
#include "../general/alnMatrixStruct.h"
#include "../general/alnSetStruct.h"
#include "../general/seqStruct.h"

/*-------------------------------------------------------\
| Fun-01 TOC: WatermanAlnNoGap
|  - Run a Waterman Smith alignment on input sequences
|    without using gap extension penalties
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
static struct alnMatrix * WatermanAlnNoGap(
    struct seqStruct *qryST, /*query sequence and data*/
    struct seqStruct *refST, /*ref sequence and data*/
    struct alnSet *settings  /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: WatermanAlnNoGap
   '  - Run a Waterman Smith alignment on input sequences
   '    without gap extension penalties
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in initial negatives for ref
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   '  o fun-01 sec-05:
   '    - Set up for returing matrix (clean up/wrap up)
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
   char *refSeqStr = refST->seqCStr + refST->offsetUL;
   char *qrySeqStr = qryST->seqCStr + qryST->offsetUL;

   /*Find the length of the reference and query*/
   ulong lenQryUL = qryST->endAlnUL - qryST->offsetUL + 1;
   ulong lenRefUL = refST->endAlnUL - refST->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   ulong lenMatrixUL = (lenRefUL + 1) * (lenQryUL + 1);
     /*+1 for the gap column and row*/

   ulong ulRef = 0;
   ulong ulQry = 0;

   /*Set up counters for the query and reference base
   `  index
   */
   /*****************************************************\
   * Fun-01 Sec-01 Sub-02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   long nextSnpScoreL = 0;/*Score for doing an match/snp*/
   long delScoreL = 0;   /*Score for doing an deletion*/

   /*Marks when to reset score buffer (every second row)*/
   long *scoreAryL = 0; /*matrix to use in alignment*/

   /*****************************************************\
   * Fun-01 Sec-01 Sub-03:
   *  - Directinol matrix variables
   \*****************************************************/

   /*Direction matrix (one cell holds a single direction)*/
   struct alnMatrix *retMatrixST = 0;
   char *dirMatrix = 0;/*Direction matrix*/
   char *insDir = 0;   /*Direction above cell*/
   ulong indexUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

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

   retMatrixST->lenRefUL = lenRefUL;
   retMatrixST->refOffsetUL = refST->offsetUL;

   retMatrixST->lenQryUL = lenQryUL;
   retMatrixST->qryOffsetUL = qryST->offsetUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Get the first indel position*/
   insDir = dirMatrix;

   for(indexUL = 0; indexUL <= lenRefUL; ++indexUL)
      dirMatrix[indexUL] = defMvStop;

   dirMatrix += indexUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Get the initial scores
   ^  o fun-01 sec-04 sub-02:
   ^    - Do the first move
   ^  o fun-01 sec-04 sub-03:
   ^    - Fill out the matrix
   ^  o fun-01 sec-04 sub-04:
   ^    - Find the next matches score
   ^  o fun-01 sec-04 sub-05:
   ^    - Find the best score for the last round
   ^  o fun-01 sec-04 sub-06:
   ^    - Find the score for the next deletion
   ^  o fun-01 sec-04 sub-07:
   ^    - Check if is an alternative base best score
   *  o fun-01 sec-04 sub-08:
   *    - Does not exist, but exists in WatermanAltAln
   ^  o fun-01 sec-04 sub-09:
   ^    - Find the scores for the next insertion
   ^  o fun-01 sec-4 sub-10:
   ^    - Move to the next reference base
   ^  o fun-01 sec-04 sub-11:
   ^    - Find the best score for the last base
   ^  o fun-01 sec-04 sub-12:
   ^    - Was the last score the best score?
   ^  o fun-01 sec-04 sub-13:
   ^    - Move to the indel column
   ^  o fun-01 sec-04 sub-14:
   ^    - Find the scores for the first base in the row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Get the initial scores
   \*****************************************************/

   nextSnpScoreL = scoreAryL[0];
   delScoreL = 0;
   dirMatrix[0] = defMvStop;

   /*These are always negative*/
   ++indexUL;
   refSeqStr = refST->seqCStr + refST->offsetUL - 1;
   qrySeqStr = qryST->seqCStr + qryST->offsetUL;

   /*****************************************************\
   * Fun-01 Sec-04 Sub-03:
   *  - Fill out the matrix
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      ulQry = 0;
      ulQry <= qryST->endAlnUL - qryST->offsetUL;
      ++ulQry
   ){ /*loop; compare query base against all ref bases*/

     for(
        ulRef = 1;
        ulRef <= refST->endAlnUL - refST->offsetUL;
        ++ulRef
     ){ /* loop; compare one query to one reference base*/

       /*************************************************\
       * Fun-01 Sec-04 Sub-04:
       *  - Find the next matches score
       \*************************************************/

       waterMaxScoreNoGap(
          refSeqStr[ulRef],
          qrySeqStr[ulQry],
          scoreAryL[ulRef],
          dirMatrix[ulRef],
          insDir[ulRef],
          nextSnpScoreL,
          delScoreL,
          settings
       );

       /*************************************************\
       * Fun-01 Sec-04 Sub-07:
       *  - Determine if is a best score (keep as primary)
       \*************************************************/

        /*This is faster than the branchless option.
        ` I am guessing to much is done in the if and 
        ` that the if if fired rarely.
        */
        if(retMatrixST->bestScoreL < scoreAryL[ulRef])
        { /*if have a new best score*/
           retMatrixST->bestScoreL = scoreAryL[ulRef];
           retMatrixST->bestEndIndexUL = indexUL;
        } /*if have a new best score*/

        ++indexUL;
     } /*loop; compare one query to one reference base*/

      /*************************************************\
      * Fun-01 Sec-04 Sub-11:
      *  - Find the best score for the last base
      \*************************************************/

      /*Find the best score for the last base. In this
      ` case, this score can only apply to indels. So,
      ` I need to move off it to avoid overwirting it
      */
      waterMaxEndRowScoreNoGap(
         refSeqStr[ulRef],
         qrySeqStr[ulQry],
         scoreAryL[ulRef],
         dirMatrix[ulRef],
         insDir[ulRef],
         nextSnpScoreL,
         delScoreL,
         settings
      );

      /***************************************************\
      * Fun-01 Sec-04 Sub-12:
      *  - Was the last score the best score?
      \***************************************************/

      /*This is faster than the branchless option.
      ` I am guessing to much is done in the if and 
      ` that the if if fired rarely.
      */
      if(retMatrixST->bestScoreL < scoreAryL[ulRef])
      { /*if have a new best score*/
         retMatrixST->bestScoreL = scoreAryL[ulRef];
         retMatrixST->bestEndIndexUL = indexUL;
      } /*if have a new best score*/

      ++indexUL;

      /**************************************************\
      *  Fun-01 Sec-04 Sub-13:
      *   - Move to the indel column
      \**************************************************/

       nextSnpScoreL = scoreAryL[0];

       insDir = dirMatrix;
       dirMatrix += ulRef + 1;
       dirMatrix[0] = defMvStop;
       scoreAryL[0] += settings->gapOpenC;
       scoreAryL[0] &= (-(scoreAryL[0] > 0));

       delScoreL = scoreAryL[0] + settings->gapOpenC;
       delScoreL &= -(delScoreL > 0);
       
       ++indexUL;
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Move back to the lower right conor cell
   ` This is not needed, but is nice.
   */

   insDir[ulRef + 1] = defMvStop;
   free(scoreAryL);
   scoreAryL = 0;

   return retMatrixST;
} /*WatermanAln*/

#endif
