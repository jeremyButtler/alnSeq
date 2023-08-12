/*#########################################################
# Name alignmentsFun
# Use:
#  o Holds functions for doing a pairwise Needleman Wunsch
#    alignment
# Libraries:
#   - "generalAlnFun.h"
#   - "alnStruct.h"
#   - "alnMatrixStruct.h"
#   o "twoBitArrays.h"
#   o "scoresST.h"
#   o "seqStruct.h"
#   o "alnSetStruct.h"
#   o "alnSeqDefaults.h"
# C Standard libraries Used:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
#########################################################*/

#include "needleman.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 NeedleManWunschAln:
'    o Perform a Needleman-Wunsch alignment on the two
'      input sequences
'  - fun-02 updateDirAndScore:
'    o Picks best score and direction for the current base
'      pairs being compared in a Needleman Wunsch alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the directoin matrix and best
|      score struct pointing to the connor left cell in
|      the direction matrix
\--------------------------------------------------------*/
struct alnMatrixStruct * NeedlemanAln(
    struct seqStruct *queryST, // query sequence and data
    struct seqStruct *refST,  // ref sequence and data
      // both queryST and refST have the sequence,
      // they also have the point to start the alignment
      // seqST->offsetUL (index 0) and the point to end
      // the alignment seqST->endAlnUL (index 0).
      // CURRENTLY THIS ONLY WORKS FOR FULL ALIGNMENTS
    struct alnSet *settings // Settings for the alignment
    // *startI and *endI paramaters should be index 1
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: NeedlemanAln
   '  - Perform a Needleman-Wunsch alignment on the two
   '    input sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Allocate memory for alignment
   '  o fun-01 sec-03:
   '    - Fill in the initial negatives for the reference
   '  o fun-01 sec-04:
   '    - Fill the matrix with scores
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Get the start of the query and reference sequences
   char *refStartCStr = refST->seqCStr + refST->offsetUL;

   char *queryStartCStr =
        queryST->seqCStr + queryST->offsetUL;

   char *tmpQueryCStr = 0;
   char *tmpRefCStr = 0;

   // Find the length of the reference and query
   unsigned long lenQueryUL =
       queryST->endAlnUL - queryST->offsetUL + 1;
       // + 1 to Account for index 0
   unsigned long lenRefUL =
       refST->endAlnUL - refST->offsetUL + 1;
       // + 1 to Account for index 0

   // Scoring variables 
   long insScoreL = 0;     // Score of the top cell
   long snpScoreL = 0; // Score of the diagnol cell
   long delScoreL = 0;    // Score of the left cell

   long *scoreMatrixL = 0; // matrix to use in alignment
   long *scoreOnLPtr = 0;  // Score I am working on
   long *lastBaseLPtr = 0; // Pointer to score above base

   // For Swaping score buffers each round
   char swapBuffBl = 1;

   // Direction matrix (one cell holds a single direction)
   struct twoBitAry *dirMatrix = 0;// Direction matrix
   struct twoBitAry topDir;        // Direction above cell
   struct twoBitAry leftDir;       // Direction before cell

   struct alnMatrixStruct *retMtxST =
     calloc(1, sizeof(struct alnMatrixStruct));

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(retMtxST == 0) return 0;

   initAlnMatrixST(retMtxST);

   scoreMatrixL =
     calloc(2 * (lenRefUL + 1) + 2, sizeof(unsigned long));
       // I need two rows to keep track of the scores (2x)
       // lenRefUL + 1 is to account for insertion column

   if(scoreMatrixL == 0)
   { // If I had a memory error
     free(retMtxST);
     return 0;
   } // If I had a memory error

   dirMatrix =
     makeTwoBitArry(
       (lenRefUL + 1) * (lenQueryUL + 1) + 2,
       0
   );
     // Calls calloc and adds an extra element at end
       // lenRefUL + 1 accounts for insertion reference row
       // lenQeurI + 1 accounts for insertion query column

   if(dirMatrix == 0)
   { // If I do not have a direction matrix for each cell
     free(retMtxST);
     free(scoreMatrixL);
     return 0;
   } // If I do not have a direction matrix for each cell

   retMtxST->dirMatrixST = dirMatrix;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in the initial negatives for the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Build up the indels for the reference row
   scoreOnLPtr = scoreMatrixL;
   *scoreOnLPtr = 0;         // Top left cell starts at 0
   ++scoreOnLPtr;            // Move to next cell (score)

   changeTwoBitElm(dirMatrix, defMvStop);
   twoBitAryMoveToNextElm(dirMatrix);

   // 2nd score (first indel in matrix)
   *scoreOnLPtr = settings->gapOpenI;
   ++scoreOnLPtr;      // Move to the 2nd reference base

   changeTwoBitElm(dirMatrix, defMvDel);
   twoBitAryMoveToNextElm(dirMatrix);

   // Set up scores for the remaning cells in the first row
   tmpRefCStr = refStartCStr + 1;

   while(*tmpRefCStr != '\0')
   { // loop; till have initalized the first row
     *scoreOnLPtr = *(scoreOnLPtr-1) +settings->gapExtendI;

     // Move to the next cell (ref base)
     ++scoreOnLPtr;
     changeTwoBitElm(dirMatrix, defMvDel);
     twoBitAryMoveToNextElm(dirMatrix);
     ++tmpRefCStr;
   } // loop; till have initalized the first row

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Set up for filling the rest of the matrix
   ^  o fun-01 sec-04 sub-02:
   ^    - Fill in the first cell (indel column)
   ^  o fun-01 sec-04 sub-03:
   ^    - Get scores for insertion, deletion, match
   ^  o fun-01 sec-04 sub-04:
   ^    - Move to the next refernce/query base
   ^  o fun-01 sec-04 sub-05:
   ^    - Handle the first cell (indel col) in new row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   ^ Fun-01 Sec-04 Sub-01:
   ^  - Set up for filling the rest of the matrix
   \******************************************************/

   // One element back will put me on the previous value
   cpTwoBitPos(dirMatrix, &leftDir);
   twoBitAryMoveBackOneElm(&leftDir);

   // The start of the limb is hte cell above the currnt
   // cell
   cpTwoBitPos(dirMatrix, &topDir);
   moveXElmFromStart(&topDir, 0);

   // Set up the previous score
   lastBaseLPtr = scoreMatrixL; // Set to first row
   swapBuffBl = 1;

   tmpQueryCStr = queryStartCStr;
   tmpRefCStr = refStartCStr;

   /******************************************************\
   * Fun-01 Sec-04 Sub-02:
   *  - Fill in the first cell (indel column)
   \******************************************************/

    *scoreOnLPtr =
      *lastBaseLPtr + settings->gapExtendI;

    changeTwoBitElm(dirMatrix, defMvIns);

    twoBitAryMoveToNextElm(dirMatrix);
    twoBitAryMoveToNextElm(&topDir);
    twoBitAryMoveToNextElm(&leftDir);

    // Move to the first base comparison
    ++scoreOnLPtr;  // Move of the indel column
    ++lastBaseLPtr; // Move off the indel column

    /*****************************************************\
    * Fun-01 Sec-04 Sub-03:
    *  - Get scores for insertion, deletion, match
    \*****************************************************/

   // Starting on the first sequence row
   while(*tmpQueryCStr != '\0')
   { // loop; fill the direction matrix with socres

       tmpRefCStr = refStartCStr;

       // Find scores & directions for all basepairs in row
       while(*tmpRefCStr != '\0')
       { // loop; compare one query to all reference bases
           snpScoreL =
              getBasePairScore(
                  tmpQueryCStr,
                  tmpRefCStr,
                  settings
           ); // Find the score for the two base pairs

           // Find the score for diagnol cell (snp/match)
           snpScoreL = *(lastBaseLPtr - 1) + snpScoreL;

           indelScore(
              insScoreL,
              getTwoBitAryElm(&topDir),
              *lastBaseLPtr,
              settings
           ); /*macro in genralAlnFun.h*/

           indelScore(
              delScoreL,
              getTwoBitAryElm(&leftDir),
              *(scoreOnLPtr - 1),
              settings
           );

           twoBitMaxScore(
             dirMatrix,       // Direction matrix
             settings,        // has direction preference
             &insScoreL,      // Score for an insertion
             &snpScoreL,  // Score for an deletion
             &delScoreL,     // Score for an match/snp
             scoreOnLPtr      // Cell on in score matrxi
           ); // Update the scores

           /**********************************************\
           * Fun-01 Sec-04 Sub-04:
           *  - Move to the next refernce/query base
           \**********************************************/
       
           // Move to the next cell to score
           ++scoreOnLPtr;  // Move to the next open score
           ++lastBaseLPtr; // Move to the score I just did
           ++tmpRefCStr;   // Move to next reference base

           // Move to the next base pair to score
           twoBitAryMoveToNextElm(dirMatrix);
           twoBitAryMoveToNextElm(&topDir);
           twoBitAryMoveToNextElm(&leftDir);
       } // loop; compare one query to all reference bases

       // Check if I still need to first row of scores
       // THis is to see what I can overwrite/reuse
       if(swapBuffBl & 1)
       { // If need to swap the buffers
           scoreOnLPtr = scoreMatrixL;
           swapBuffBl = 0;
       } // If need to swap the buffers

       else // Else 2nd row of scores is no longer needed
       { // Else need to reset the score part
           lastBaseLPtr = scoreMatrixL;
           swapBuffBl = 1;
       } // Else need to reset the score part

       ++tmpQueryCStr; // Move to the next query base

       /**************************************************\
       *  Fun-01 Sec-04 Sub-05:
       *    - Handle the first cell (indel col) in new row
       \**************************************************/

       *scoreOnLPtr =
         *lastBaseLPtr + settings->gapExtendI;

       changeTwoBitElm(dirMatrix, defMvIns);
       twoBitAryMoveToNextElm(dirMatrix);
       twoBitAryMoveToNextElm(&topDir);
       twoBitAryMoveToNextElm(&leftDir);

       // Move to the first base comparison
       ++scoreOnLPtr;  // Move to first base pair in row
       ++lastBaseLPtr; // Move to last base pair
   } // loop; fill the direction matrix with socres

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Account for being two cells out of bounds
   twoBitAryMoveBackOneElm(dirMatrix);
   changeTwoBitElm(dirMatrix, defMvStop);
   twoBitAryMoveBackOneElm(dirMatrix);

   // Set the best score to the conor right cell
   if(swapBuffBl != 0)
     retMtxST->bestScoreST.scoreL = *(scoreOnLPtr - 2);
   else
     retMtxST->bestScoreST.scoreL =
       *(lastBaseLPtr + lenRefUL - 1);

   retMtxST->bestScoreST.indexUL =
     getTwoBitAryIndex(dirMatrix);

   free(scoreMatrixL);

   return retMtxST;
} // NeeldeManWunschAln
