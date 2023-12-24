/*########################################################
# Name memWater
# Use:
#  o Holds functions doing a memory efficent Smith
#    Waterman pairwise alignments. These aligners return
#    indexes, which can be convereted and then run through
#    used to find a  global/local aligment.
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
#   o <stdio.h>
########################################################*/

#ifndef MEM_WATER_H
#define MEM_WATER_H

#include "../general/genScan.h"
#include "../general/alnSetStruct.h"
#include "../general/seqStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memWater SOF: Start Of Functions
' o fun-01 memWater:
'   - Run a memory efficent Waterman Smith alignment on
'     input sequences
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Fun-01: memWater
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
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
|    o alnMatrix struct with the best score
|    o 0 for memory allocation errors
\-------------------------------------------------------*/
static struct alnMatrix * memWater(
    struct seqStruct *qryST, /*query sequence*/
    struct seqStruct *refST, /*ref sequence*/
    struct alnSet *settings  /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: memWaterAln
   '  - Run a memory efficent Waterman Smith alignment on
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

   /*Iterators for loops*/
   ulong ulRefBase = 0;
   ulong ulQryBase = 0;

   /*****************************************************\
   * Fun-01 Sec-01 Sub-02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   long delScoreL = 0;    /*Score for doing an deletion*/
   long nextSnpScoreL = 0;/*Score for the next match/snp*/
   long *scoreAryL = 0;  /*matrix to use in alignment*/

   /*Used in finding if useing gap extension*/
   short gapDiffS =
      settings->gapExtendC - settings->gapOpenC;

   /*****************************************************\
   * Fun-01 Sec-01 Sub-03:
   *  - Directional matrix variables
   \*****************************************************/

   /*Direction matrix (1 cell holds a single direction)*/
   struct alnMatrix *retMatrixST = 0;
   char *dirRow = 0;  /*Holds directions*/

   /*Keeping track of alignment starting positions*/
   ulong indexUL = 0;     /*Index I am at in the matrix*/
   ulong *indexAryUL=0;   /*Row of starting indexes*/
   ulong *oldIndexAryUL=0;/*Last round starting indexes*/
   ulong *swapPtrUL = 0;  /*For swapping ulongs*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   ^  o fun-01 sec-02 sub-01:
   ^    - Allocate memory for the alignment
   ^  o fun-01 sec-02 sub-02:
   ^    - Allocate memory for keeping track of indexes
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Allocate memory for the alignment
   \****************************************************/

   retMatrixST = malloc(sizeof(struct alnMatrix));
   if(retMatrixST == 0) return 0;
   initAlnMatrix(retMatrixST);

   dirRow = malloc((lenRefUL + 1) * sizeof(char));

   if(dirRow == 0)
   { /*If: Memory error*/
      freeAlnMatrix(retMatrixST);
      return 0;
   } /*If: Memory error*/

   scoreAryL = malloc((lenRefUL + 1) * sizeof(long));
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
      indexAryUL[indexUL] = indexUL;
      scoreAryL[indexUL] = 0;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Final set up before scoring the matrix
   ^  o fun-01 sec-04 sub-02:
   ^    - Start loops and get each score
   ^  o fun-01 sec-04 sub-03:
   ^    - Check if is an alternative base best score
   ^  o fun-01 sec-04 sub-04:
   ^    - Find the best score for the last base
   ^  o fun-01 sec-04 sub-05:
   ^    - Is last base in row an alternative alignment?
   ^  o fun-01 sec-04 sub-06:
   ^    - Prepare to score the next row in the matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Final set up before scoring the matrix
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
   *  - Start loops and get each score
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      ulQryBase = 0;
      ulQryBase <= qryST->endAlnUL - qryST->offsetUL;
      ++ulQryBase
   ){ /*loop; compare query base against all ref bases*/

      for(
         ulRefBase = 1;
         ulRefBase < lenRefUL;
         ++ulRefBase
      ){ /*loop; compare one query to one reference base*/
         waterScanMaxScore(
            refSeqStr[ulRefBase],
            qrySeqStr[ulQryBase],
            gapDiffS,
            scoreAryL[ulRefBase],
            dirRow[ulRefBase],
            indexAryUL[ulRefBase],
            nextSnpScoreL,
            delScoreL,
            oldIndexAryUL[ulRefBase - 1], 
            oldIndexAryUL[ulRefBase], 
            indexAryUL[ulRefBase - 1],
            indexUL,
            settings
         );

         /***********************************************\
         * Fun-01 Sec-04 Sub-03:
         *  - Determine if is best score (keep as primary)
         \***********************************************/

         if(retMatrixST->bestScoreL <scoreAryL[ulRefBase])
         { /*If: this is the best score*/
            retMatrixST->bestScoreL =scoreAryL[ulRefBase];

            retMatrixST->bestStartIndexUL =
               indexAryUL[ulRefBase];

            retMatrixST->bestEndIndexUL = indexUL;
         } /*If: this was an snp or match*/

         ++indexUL;
      } /*loop; compare one query to one reference base*/

      /***************************************************\
      * Fun-01 Sec-04 Sub-04:
      *  - Find the best score for the last base
      \***************************************************/

      waterScanMaxEndRowScore(
         refSeqStr[ulRefBase],
         qrySeqStr[ulQryBase],
         gapDiffS,
         scoreAryL[ulRefBase],
         dirRow[ulRefBase],
         indexAryUL[ulRefBase],
         nextSnpScoreL,
         delScoreL,
         oldIndexAryUL[ulRefBase - 1], 
         oldIndexAryUL[ulRefBase], 
         indexAryUL[ulRefBase - 1],
         indexUL,
         settings
      );

     /***************************************************\
     * Fun-01 Sec-04 Sub-05:
     *  - Is last base in row an alternative alignment?
     \***************************************************/

     /*This is one part were a branched operation is
     ' more efficent. I think this is because this
     ' branch only invovles a quick if check
     */
     if(retMatrixST->bestScoreL < scoreAryL[ulRefBase])
     { /*If: this is the best score*/
        retMatrixST->bestScoreL = scoreAryL[ulRefBase];

        retMatrixST->bestStartIndexUL =
           indexAryUL[ulRefBase];

        retMatrixST->bestEndIndexUL = indexUL;
     } /*If: this was an snp or match*/

     ++indexUL;

     /***************************************************\
     *  Fun-01 Sec-04 Sub-06:
     *   - Prepare for the next round
     \***************************************************/

      /*Get scores set up*/
      nextSnpScoreL = scoreAryL[0];

      dirRow[0] = defMvStop;
      scoreAryL[0] += settings->gapExtendC;
      scoreAryL[0] &= -(scoreAryL[0] > 0);

      delScoreL = scoreAryL[0] + settings->gapExtendC;
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
   free(scoreAryL);
   free(indexAryUL);
   free(oldIndexAryUL);

   dirRow = 0;
   scoreAryL = 0;
   indexAryUL = 0;
   oldIndexAryUL = 0;

   return retMatrixST;
} /*memWaterAln*/

#endif
