/*########################################################
# Name memWaterScanNoGap
# Use:
#  o Holds functions doing a memory efficent Smith
#    Waterman pairwise alignments with a query reference
#    scan. These functions return a set of indexes, which
#    can be convereted and then run through used to find a
#    global/local aligment.
#  o This variation trys to get a best score for each
#    query and reference base and does not use the
#    gap extension penalties
#  o Coversion macros are in alnMatrixsStruct.h:
#    indexToQry (fun-07), indexToRef (fun-08), and
#    indexToCoord (fun-09)
# Libraries:
#   - "../general/genScan.h"           (No .c file)
#   o "../general/alnMatrixStruct.h"   (No .c file)
#   o "../general/twoBitArrays.h"      (No .c file)
#   o "../general/dataTypeShortHand.h" (No .c file)
#   o "../general/genMath.h"           (No .c file)
#   - "../general/alnSetStruct.h"      (No .c file)
#   o "../general/alnSeqDefaults.h"    (No .c file)
#   - "../general/seqStruct.h"         (No .c file)
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
########################################################*/

#ifndef MEM_WATER_SCAN_NO_GAP_H
#define MEM_WATER_SCAN_NO_GAP_H

#include "../general/genScan.h"
#include "../general/alnSetStruct.h"
#include "../general/seqStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memWaterAltAln SOF: Start Of Functions
' o fun-01 memWaterAltAln:
'   - Run a memory efficent Waterman Smith alignment that
'     returns alternative alignmetns
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Fun-01: memWaterScanNoGap
|   - Performs a memory efficent Smith Waterman alignment
|     with a query/reference scan on a pair of sequences
| Input;
|   - qryST:
|     o SeqStruct with the query sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment.
|   - refST:
|     o SeqStruct with the reference sequence and index 0
|       coordinates to start (offsetUL)/end (endAlnUL) the
|       alignment
|   - settings:
|     o alnSet structure with the setttings, such as
|       gap open, gap extend, scoring matrix, and
|       preffered direction.
| Output:
|  - Returns:
|    o an alnMatrix struct with the best scores
|    o 0 for memory allocation errors
\-------------------------------------------------------*/
static struct alnMatrix * memWaterScanNoGap(
    struct seqStruct *qryST, /*query sequence*/
    struct seqStruct *refST, /*ref sequence*/
    struct alnSet *settings  /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: memWaterScanNoGap
   '  - Run a memory efficent Waterman Smith alignment
   '    scan without gap extension penalties on the
   '    input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in initial negatives for ref
   '  o fun0 sec-04:
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
   char *refSeqStr = 0;
   char *qrySeqStr = 0;

   ulong lenRefUL = refST->endAlnUL - refST->offsetUL + 1;
   ulong lenQryUL = qryST->endAlnUL - qryST->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   /*Iterators for loops*/
   ulong ulRef = 0;
   ulong ulQry = 0;

   /*****************************************************\
   * Fun-01 Sec-01 Sub-02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   long delScoreL = 0;    /*Score for doing an deletion*/
   long nextSnpScoreL = 0;/*Score for the next match/snp*/

   /*Marks when to reset score buffer (every second row)*/
   long *scoreAryL = 0; /*matrix to use in alignment*/

   /*****************************************************\
   * Fun-01 Sec-01 Sub-03:
   *  - Directional matrix variables
   \*****************************************************/

   struct alnMatrix *retMatrixST = 0;
   char *dirRow = 0;  /*Holds directions*/

   /*For recording alternative alignments*/
   ulong indexUL = 0;      /*Index I am at in the matrix*/
   ulong *indexAryUL=0;    /*Row of starting indexes*/
   ulong *oldIndexAryUL=0;/*Last round starting indexes*/
   ulong *swapPtrUL = 0;      /*For swapping ulongs*/

   /*For keeping track of best scores*/
   ulong lenAltRowUL = 0;

   long *refScoreAryL = 0;
   ulong *refIndexAryUL = 0;
   ulong *refEndIndexAryUL = 0;

   long *qryScoreAryL = 0;
   ulong *qryIndexAryUL = 0;
   ulong *qryEndIndexAryUL = 0;

   ulong ulScore = 0; /*For for loop at end*/

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

   retMatrixST = malloc(1 * sizeof(struct alnMatrix));
   if(retMatrixST == 0) return 0;
   initAlnMatrix(retMatrixST);

   dirRow = malloc((lenRefUL + 1) * sizeof(char));

   if(dirRow == 0)
   { /*If: Memory error*/
      freeAlnMatrix(retMatrixST);
      free(dirRow);
      return 0;
   } /*If: Memory error*/

   scoreAryL = calloc((lenRefUL + 1), sizeof(long));
   /*+ 1 is for the indel column*/

   if(scoreAryL == 0)
   { /*If I had a memory error*/
     freeAlnMatrix(retMatrixST);
     free(dirRow);
     return 0;
   } /*If I had a memory error*/

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
   refScoreAryL = calloc(lenAltRowUL, sizeof(long));

   if(refScoreAryL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      free(dirRow);
      return 0;
   } /*If had a memory error*/

   /*I am running the loops with the
   ` query at index 0 and the reference
   ` at index 1. So I need a + 1
   */
   retMatrixST->scoreAryL = refScoreAryL;
   qryScoreAryL = refScoreAryL + lenRefUL + 1;

   /*Set up the row of starting indexes*/
   refIndexAryUL = calloc(lenAltRowUL, sizeof(ulong));

   if(refIndexAryUL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      free(dirRow);
      return 0;
   } /*If had a memory error*/

   retMatrixST->startIndexAryUL = refIndexAryUL;
   qryIndexAryUL = refIndexAryUL + lenRefUL + 1;

   /*Set up the row of ending indexes*/
   refEndIndexAryUL = calloc(lenAltRowUL, sizeof(ulong));

   if(refEndIndexAryUL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      free(dirRow);
      return 0;
   } /*If had a memory error*/

   retMatrixST->endIndexAryUL = refEndIndexAryUL;
   qryEndIndexAryUL = refEndIndexAryUL + lenRefUL + 1;

   /*****************************************************\
   * Fun-05 Sec-02 Sub-03:
   *  - Get memory for keeping track of starting indexes
   \*****************************************************/

   /*Set up the first row of starting indexes*/
   indexAryUL = malloc((lenRefUL + 1) * sizeof(ulong));

   if(indexAryUL == 0)
   { /*If had a memory error*/
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      free(dirRow);
      return 0;
   } /*If had a memory error*/

   /*Set up the second row of indexs (so have two rows)*/
   oldIndexAryUL = malloc((lenRefUL + 1) * sizeof(ulong));

   if(oldIndexAryUL == 0)
   { /*If had a memory error*/
      free(indexAryUL);
      free(scoreAryL);
      freeAlnMatrix(retMatrixST);
      free(dirRow);
      return 0;
   } /*If had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(indexUL = 0; indexUL <= lenRefUL; ++indexUL)
   { /*loop; till have initalized the first row*/
      dirRow[indexUL] = defMvStop;
      indexAryUL[indexUL] = indexUL; /*Stop is new index*/
   } /*loop; till have initalized the first row*/
   /*Everthing else was initalized by calloc*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Final setup for the scoring loop
   ^  o fun-01 sec-04 sub-02:
   ^    - Start loops and get scores for the first half
   ^      of the row (reference bases)
   ^  o fun-01 sec-04 sub-03:
   ^    - Check if keeping score for first half of row
   ^  o fun-01 sec-04 sub-04:
   ^    - Find the scores for the second half of each row
   ^  o fun-01 sec-04 sub-05:
   ^    - Check if kedeping score for 2nd half of each row
   ^  o fun-01 sec-04 sub-06:
   ^    - Find the score for the last base in each row
   ^  o fun-01 sec-04 sub-07:
   ^    - Check if keeping the last score in each row
   ^  o fun-01 sec-04 sub-08:
   ^    - Setup for the next loop
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Finall set up for the scoring loop
   \*****************************************************/

   /*Move the row of starting indexes to the last row*/
   swapPtrUL = indexAryUL;
   indexAryUL = oldIndexAryUL;
   oldIndexAryUL = swapPtrUL;

   nextSnpScoreL = scoreAryL[0];

   /*These are always negative*/
   delScoreL = 0;
   indexAryUL[0] = indexUL;
   dirRow[0] = defMvStop;

   /*Incurment to the frist base*/
   ++indexUL;
   refSeqStr = refST->seqCStr + refST->offsetUL - 1;
   qrySeqStr = qryST->seqCStr + qryST->offsetUL;

   /*****************************************************\
   * Fun-01 Sec-04 Sub-02:
   *  - Start loops and get each score for the frist half
   *    of each row
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      ulQry = 0;
      ulQry <= qryST->endAlnUL - qryST->offsetUL;
      ++ulQry
   ){ /*loop; compare query base against all qry bases*/

      for(
         ulRef = 1;
         ulRef < (lenRefUL)/ 2;
         ++ulRef
      ){ /*loop; compare one query to one reference base*/
         waterScanMaxScoreNoGap(
            refSeqStr[ulRef],
            qrySeqStr[ulQry],
            scoreAryL[ulRef],
            dirRow[ulRef],
            indexAryUL[ulRef],
            nextSnpScoreL,
            delScoreL,
            oldIndexAryUL[ulRef - 1],
            oldIndexAryUL[ulRef],
            indexAryUL[ulRef - 1],
            indexUL,
            settings
         );

         /***********************************************\
         * Fun-01 Sec-04 Sub-03:
         *  - Check if keeping score for first half of
         *    row
         \***********************************************/

         scanIfKeepScoreRef(
            scoreAryL[ulRef],
            dirRow[ulRef],
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
      } /*loop; compare one query to one reference base*/

      /***********************************************\
      * Fun-01 Sec-04 Sub-04:
      *  - Find scores for the second half of each row
      \***********************************************/

      /*I am doing this in two loops so I can have the
      ` reference get priotority for scores in the first
      ` quadrent and the query gets priority for scores
      ` in the upper right quadrent
      */

      for(
         ulRef = ulRef;
         ulRef < lenRefUL;
         ++ulRef
      ){ /*loop; compare one query to one reference base*/
         waterScanMaxScoreNoGap(
            refSeqStr[ulRef],
            qrySeqStr[ulQry],
            scoreAryL[ulRef],
            dirRow[ulRef],
            indexAryUL[ulRef],
            nextSnpScoreL,
            delScoreL,
            oldIndexAryUL[ulRef - 1],
            oldIndexAryUL[ulRef],
            indexAryUL[ulRef - 1],
            indexUL,
            settings
         );

         /***********************************************\
         * Fun-01 Sec-04 Sub-05:
         *  - Check if keeping score for 2nd half of row
         \***********************************************/

         scanIfKeepScoreQry(
            scoreAryL[ulRef],
            dirRow[ulRef],
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
      } /*loop; compare one query to one reference base*/

      /***************************************************\
      * Fun-01 Sec-04 Sub-06:
      *  - Find the best score for last base of each row
      \***************************************************/

      waterScanMaxEndRowScoreNoGap(
         refSeqStr[ulRef],
         qrySeqStr[ulQry],
         scoreAryL[ulRef],
         dirRow[ulRef],
         indexAryUL[ulRef],
         nextSnpScoreL,
         delScoreL,
         oldIndexAryUL[ulRef - 1],
         oldIndexAryUL[ulRef],
         indexAryUL[ulRef - 1],
         indexUL,
         settings
      );

     /***************************************************\
     * Fun-01 Sec-04 Sub-07:
     *  - Is last base in row an alternative alignment?
     \***************************************************/

     scanIfKeepScoreQry(
        scoreAryL[ulRef],
        dirRow[ulRef],
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
     *  Fun-01 Sec-04 Sub-08:
     *   - Prepare for the next round
     \***************************************************/

      /*Get scores set up*/
      nextSnpScoreL = scoreAryL[0];
      dirRow[0] = defMvStop;
      scoreAryL[0] += settings->gapOpenC;
      scoreAryL[0] &= -(scoreAryL[0] > 0);

      delScoreL = scoreAryL[0] + settings->gapOpenC;
      delScoreL &= -(delScoreL > 0);

      /*Swap index arrays so the current is last*/
      swapPtrUL = indexAryUL;
      indexAryUL = oldIndexAryUL;
      oldIndexAryUL = swapPtrUL;

      ++indexUL; /*Set index for the next base pair*/
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

   free(dirRow);
   free(indexAryUL);
   free(oldIndexAryUL);
   free(scoreAryL);

   dirRow = 0;
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
} /*memWaterScanNoGap*/

#endif
