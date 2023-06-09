/*#########################################################
# Name alignmentsFun
# Use:
#  o Holds functions for doing pairwise alignments (Needleman/waterman)
# Includes:
#   - "defaultSettings.h"
#   - "cStrToNumberFun.h"
#   - "twoBitArrays.h"
# C Standard libraries:
#   - <stdlib.h>
#   - <stdio.h>
#   o <stdint.h>
#########################################################*/

#include "needleman.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 NeedleManWunschAln:
'     o Perform a Needleman-Wunsch alignment on input sequences
'  - fun-02 updateDirAndScore:
'     o Picks the best score and direction for the current base pairs
'       being compared in a Needleman Wunsch alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*---------------------------------------------------------------------\
| Output:
|  - Returns:
|    o array with flags for snp/match, insertion, and deletions at each
|      position. (1 = snp/match, 2 = insertion, 4 = deletion)
|    o 0 for memory allocation errors
|  - Modifies:
|    o lenErrAryUI to hold the length of the returned array
\---------------------------------------------------------------------*/
struct alnMatrixStruct * NeedlemanAln(
    struct seqStruct *queryST,
      // query sequence, length, & bounds for alignment
    struct seqStruct *refST,  
      // reference sequence, length, & bounds for alignment
    struct alnSet *settings // Settings for the alignment
    // *startI and *endI paramaters should be index 1
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: NeedlemanAln
   '  - Perform a Needleman-Wunsch alignment on input sequences
   '  o fun-01 sec-1: Variable declerations
   '  o fun-01 sec-2: Allocate memory for alignment
   '  o fun-01 sec-3: Fill in the initial negatives for the reference
   '  o fun-01 sec-4: Fill the matrix with scores
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-1: Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char firstRoundBl = 1;
   char swapBuffBl = 1;   // Swap buffers when finshed every 2nd row

   // Get the start of the query and reference sequences
   char *refStartCStr =
       refST->seqCStr
     + refST->offsetUI
     - 1;

   char *queryStartCStr =
        queryST->seqCStr
      + queryST->offsetUI
      - 1;

   char *tmpQueryCStr = queryStartCStr;
   char *tmpRefCStr = refStartCStr;

   // Find the length of the reference and query
   unsigned long lenQueryUL =
       queryST->endAlnUI
     - queryST->offsetUI
     + 1; // +1 converts to index 1 (subtraction makes 0)

   unsigned long lenRefUL =
       refST->endAlnUI
     - refST->offsetUI
     + 1; // +1 converts to index 1 (subtraction makes 0)

   unsigned long lenMatrixUL = 0;

       // +1 to convert back to 1 index (subtraction makes 0)

   long snpScoreL = 0;     // Score for single base pair
   long scoreTopL = 0;     // Score when using the top cell
   long scoreDiagnolL = 0; // Score when using the diagnol cell
   long scoreLeftL = 0;    // Score when using the left cell

   long *scoreMatrixL = 0; // matrix to use in alignment
   long *scoreOnLPtr = 0;  // Score I am currently working on
   long *lastBaseLPtr = 0; // Pointer to cell with last base

   struct twoBitAry *dirMatrixUC = 0;  // Directions for each score cell
   uint8_t *dirOnUCPtr = 0;   // Score working on

   struct twoBitAry topDir;
   struct twoBitAry leftDir;   // Deletion direction

   struct alnMatrixStruct *retMtxST =
     calloc(1, sizeof(struct alnMatrixStruct));

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-2: Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(retMtxST == 0) return 0;

   initAlnMatrixST(retMtxST);

   scoreMatrixL =
     calloc(2 * (lenRefUL + 1), sizeof(unsigned long));
       // I need two rows to keep track of the scores (2x)
       // lenRefUL + 1 is to account for insertion zero query column

   if(scoreMatrixL == 0)
   { // If I had a memory error
     free(retMtxST);
     return 0;
   } // If I had a memory error

   dirMatrix =
     makeTwoBitArry((lenRefUL + 1) * (lenQueryUL + 1), 0);
     // Calls calloc and adds an extra element at end
       // lenRefUL + 1 accounts for insertion reference row
       // lenQeurI + 1 accounts for insertion query column


   if(dirMatrixUC == 0)
   { // If I do not have a direction matrix for each cell
     free(retMtxST);
     free(scoreMatrixL);
     return 0;
   } // If I do not have a direction matrix for each cell

   retMtxST->dirMatrixST = dirMatrix;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-3: Fill in the initial negatives for the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Build up the indels for the reference row
   scoreOnLPtr = scoreMatrixL;
   *scoreOnLPtr = 0;         // Top left cell
   ++scoreOnLPtr;            // Move to first reference base cell
   *scoreOnLPtr = settings->gapStartPenaltyI;
       // Second column of the first row holds the first indel
   ++scoreOnLPtr;            // Move to first reference base cell

   // Get direction matrix in sync with scoring matrix
   changeTwoBitElm(dirMatrix, defMoveStop);
   twoBitAryMoveToNextElm(dirMatrix);

   changeTwoBitElm(dirMatrix, defMoveLeft);
   twoBitAryMoveToNextElm(dirMatrix);

   tmpRefCStr = refStartCStr + 1; // Move off the first base (is done)

   // Already filled in two cells in this row so, it is lenRef - 1
   while(*tmpRefCStr != '\0')
   { // loop; till have initalized the first row
       *scoreOnLPtr =
         *(scoreOnLPtr - 1) + settings->gapExtendPenaltyI;

       ++scoreOnLPtr; // Move to the next element
       // Move past elements in the bit array
       // Not worried about specifiying left shift (0),
       // since calloc already set everything to 0.

       changeTwoBitElm(dirMatrix, defMoveLeft);
       twoBitAryMoveToNextElm(dirMatrix);
       ++tmpRefCStr;
   } // loop; till have initalized the first row

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-4: Fill the matrix with scores
   ^  o fun-01 sec-4 sub-1: Fill in the indel column
   ^  o fun-01 sec-4 sub-2: Get scores for insertion, deletion, match
   ^  o fun-01 sec-4 sub-3: Move to the next refernce/query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Fun-01 Sec-4 Sub-1: Fill in the indel column
   \*******************************************************************/

   // Marks the last cell in the score matrix (just to make sure)
   cpTwoBitPos(dirMatrix, &leftDir);
   twoBitAryMoveBackOneElm(&leftDir);
   cpTwoBitPos(dirMatrix, &topDir);

   firstRoundBl = 1; // Mark need to add a gap start penalty at start
   swapBuffBl = 1;   // Swap buffers when finsh every 2nd row

   tmpQueryCStr = queryStartCStr;
   tmpRefCStr = refStartCStr;

   // Starting on the first sequence row
   while(*tmpQueryCStr != '\0')
   { // loop; compare one query base against all reference bases

       switch(firstRoundBl)
       { // Switch, check if on the 1st cell of the 2nd row
           case 1:
             *scoreOnLPtr =
               *lastBaseLPtr + settings->gapStartPenaltyI;
             firstRoundBl = 0; // Never will be reset
             break;

           case 0:
             *scoreOnLPtr =
               *lastBaseLPtr + settings->gapExtendPenaltyI;
             break;
           // Else is the first indel for the indel column
       } // Switch, check if on the 1st cell of the 2nd row
 
       changeTwoBitElm(dirMatrix, defMoveUp);
       twoBitAryMoveToNextElm(dirMatrix);
       twoBitAryMoveToNextElm(&topDir);
       twoBitAryMoveToNextElm(&leftDir);

       // Move to the first base comparison
       ++scoreOnLPtr;  // Get of negative column for the new query base
       ++lastBaseLPtr; // Get of negative column for last query base

       tmpRefCStr = refStartCStr;

       /***************************************************************\
       * Fun-01 Sec-4 Sub-2: Get scores for insertion, deletion, match
       \***************************************************************/

       // First reference bases column (indel column already handled)
       while(*tmpRefCStr != '\0')
       { // loop; compare one query to one reference base

           snpScoreL =
              getBasePairScore(
                  tmpQueryCStr,
                  tmpRefCStr,
                  settings
           ); // Find the score for the two base pairs

           // Find the score for the diagnol cell (snp/match)
           scoreDiagnolL = *(lastBaseLPtr - 1) + snpScoreL;

           switch(settings->matchPriorityBl) 
           { // Switch: check if matches have priority

               case 1:
                   if(
                     checkIfBasesMatch(
                       tmpQueryCStr,
                       tmpRefCStr
                     ) != 0
                   ) { // If had matching bases
                       changeTwoBitElm(
                         dirMatrix,
                         defMoveDiagnol
                       );
                       *scoreOnLPtr = scoreDiagnolL;
                       break;
                   } // If had matching bases

               case 0:   // Either not using match priority or not match
                   scoreTopL =
                       getIndelScore(
                           topDirUCPtr,
                           &topBitUC,
                           settings,
                           lastBaseLPtr
                   ); // Get the score for an insertion

                   // If the limb is not complete the last direction
                   // will always be one shift back 
                   if(leftBitUC < 3) tmpLeftBitUC = 2;
                   else tmpLeftBitUC = 3;

                   scoreLeftL =
                       getIndelScore(
                           leftDirUCPtr,   // part of two bit index
                           &tmpLeftBitUC,  // part of two bit index
                           settings,       // Has gap penalties
                           scoreOnLPtr - 1 // Score of the previous base
                   ); // Get the score for an insertion

                   updateDirAndScoreNeedle(
                       dirOnUCPtr,
                       settings,   // Has preference for score selection
                       &scoreTopL,     // Score for an insertion
                       &scoreDiagnolL, // Score for an match/snp
                       &scoreLeftL,    // The score for an deletion
                       scoreOnLPtr
                   ); // Update the scores
           } // Switch: check if matches have priority

           /************************************************************\
           * Fun-01 Sec-4 Sub-3: Move to the next refernce/query base
           \************************************************************/
       
           // Move to the next cell to score
           ++scoreOnLPtr; // Move to next comparison for this query base
           ++lastBaseLPtr; // Move to next element
           ++tmpRefCStr;   // Move to the next reference base

           // Move to the next base pair to score
           twoBitAryMoveToNextElm(dirMatrix);
           twoBitAryMoveToNextElm(topDir);
           twoBitAryMoveToNextElm(leftDir);
       } // loop; compare one query to one reference base

       // Will end on the corrnor cell
       *scoreL = *(scoreOnLPtr - 1);

       if(swapBuffBl & 1)
       { // If need to swap the buffers
           scoreOnLPtr = scoreMatrixL; // Restarting scoring
           swapBuffBl = 0;
       } // If need to swap the buffers

       else
       { // Else need to reset the score part
           lastBaseLPtr = scoreMatrixL; // Last base is on first row
           swapBuffBl = 1;
       } // Else need to reset the score part

       ++tmpQueryCStr; // Move to the next query base
   } // loop; compare one query base against all reference bases

   return retMtxST;
} // NeeldeManWunschAln

/*---------------------------------------------------------------------\
| Output: Modifies: scoreOnL and dirOnUC to hold best score & direction
\---------------------------------------------------------------------*/
void updateDirAndScoreNeedle(
    uint8_t *dirOnUCPtr,     // Direction on with first two bits cleared
    struct alnSet *alnSetST, // Has preference for score selection
    long *scoreTopL,     // Score for an insertion
    long *scoreDiagnolL, // Score for an match/snp
    long *scoreLeftL,    // The score for an deletion
    long *scoreOnL       // Score to update
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: updateDirAndScoreNeedle
   '  - Picks the best score and direction for the current base pairs
   '    being compared in a Needleman Wunsch alignment
   '  o fun-02 sec-1: Matches->insertions->deletions
   '  o fun-02 sec-2: Matches->deletions->insertions
   '  o fun-02 sec-3: Insertions->matches->deletions
   '  o fun-02 sec-4: Deletions->matches->insertions
   '  o fun-02 sec-5: Insertions->deletions->matches
   '  o fun-02 sec-6: Deletions->insertions->matches
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // I decide the best path as I score since all equal scores create
   // an equally valid alignment.

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-1: Matches and then insertions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   switch(alnSetST->diagnolPriorityC)
   { // Switch; get an snp/match priority
       case 0:  // Top priority
       // Case; bases or matches are top priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 1:
               // Case: priority is matches/snps and then insertions
                   if(*scoreDiagnolL >= *scoreTopL)
                   { // If diagnol beats insertion

                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats deletions

                       else
                       { // Else the deletion is the best score
                           *dirOnUCPtr |= defMoveLeft;
                           *scoreOnL = *scoreLeftL;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score

                   else
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score
                       
                   return;
               // Case: priority is matches/snps and then insertions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-02 Sec-2: Matches->deletions->insertions
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is matches/snps and then deletions
                   if(*scoreDiagnolL >= *scoreLeftL)
                   { // If diagnol beats deletion

                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the insertion is the best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score
                       
                   return;
               // Case: priority is matches/snps and then deletions
           } // Priority for insertions
       // Case; bases or matches are top priority

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-02 Sec-3: Insertions->matches->deletions
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 1:
       // Case; bases or matches are second priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority is insertions and then matches/snps
                   if(*scoreDiagnolL > *scoreTopL)
                   { // If diagnol beats insertion

                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats deletions

                       else
                       { // Else the deletion is the best score
                           *dirOnUCPtr |= defMoveLeft;
                           *scoreOnL = *scoreLeftL;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score

                   else
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score
                       
                   return;
               // Case: priority is matches/snps and then insertions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-02 Sec-4: Deletions->matches->insertions
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is deletions and then matches/snps
                   if(*scoreDiagnolL > *scoreLeftL)
                   { // If diagnol beats deletion

                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the insertion is the best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score
                       
                   return;
               // Case: priority is matches/snps and then deletions
           } // Priority for insertions

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-02 Sec-5: Insertions->deletions->matches
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 2:
       // Case; bases or matches are last priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority is insertions and then deletions
                   if(*scoreTopL >= *scoreLeftL)
                   { // If diagnol beats insertion

                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol is the highest score
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol is the highest score

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the match/snp is the best score
                      *dirOnUCPtr |= defMoveDiagnol;
                      *scoreOnL = *scoreDiagnolL;
                   } // Else the match/snp is the best score
                       
                   return;
               // Case: priority is insertions then deletions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-02 Sec-6: Deletions->insertions->matches
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is deletions and then insertions
                   if(*scoreTopL > *scoreLeftL)
                   { // If insertion beats deletion

                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol beats insertions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the insertion is the best score
                   } // If insertion beats deletion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the match/snp is the best score
                      *dirOnUCPtr |= defMoveDiagnol;
                      *scoreOnL = *scoreDiagnolL;
                   } // Else the match/snp is the best score
                       
                   return;
               // Case: priority is insertions and then deletions
           } // Priority for insertions
       // Case; bases or matches are last priority
   } // Switch; get an snp/match priority

   return;
} // updateDirAndScoreNeedle

