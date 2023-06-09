/*#########################################################
# Name waterman
# Use:
#  o Holds functions doing a Waterman-Smith pairwise
#    alignments. This version outputs a single alignment
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
#   o <stdio.h>  // by alnSetStructure.h
#   - <string.h>
#########################################################*/

#include "waterman.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 WatermanSmithAln:
'     o Perform a Waterman Smith alignment on input
'       sequences
'  o fun-02 addBestBaseScore:
'    - Adds a score and index to the kept scores list
'  o fun-03 printMatrixCig:
'    - Prints out a cigar for an single path in a
'      direction matrix
'  o fun-04 printAltWaterAlns:
'    - Prints out the best aligment and the saved
'       alterantive alignments  (best alignment for each
'       base) to a file
'  - fun-05 updateDirScoreWaterSingle:
'     o Picks the best score and direction for the current
'       base pairs being compared in a Waterman-Smith
'       alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnMatrixStruct with the direction matrix and scores
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
struct alnMatrixStruct * WatermanAln(
    struct seqStruct *queryST, // query sequence and data
    struct seqStruct *refST,  // ref sequence and data
      // both queryST and refST have the sequence,
      // they also have the point to start the alignment
      // seqST->offsetUL (index 0) and the point to end
      // the alignment seqST->endAlnUL (index 0).
      // CURRENTLY THIS ONLY WORKS FOR FULL ALIGNMENTS
    struct alnSet *setST,// Settings for the alignment
    // *startI and *endI paramaters should be index 1
    char *prefixCStr  // Prefix for matrix scan output file
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: WatermanAln
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
   '    - Set up for returing the matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
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
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-01 Sub-01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \******************************************************/

   // Get the start of the query and reference sequences
   char *refStartCStr = refST->seqCStr + refST->offsetUL;

   char *queryStartCStr =
        queryST->seqCStr + queryST->offsetUL;

   char *tmpQueryCStr = queryStartCStr;
   char *tmpRefCStr = refStartCStr;

   // Find the length of the reference and query
   unsigned long lenQueryUL =
       queryST->endAlnUL - queryST->offsetUL + 1;
     // The + 1 is to account for index 0 of endAlnUL

   unsigned long lenRefUL =
       refST->endAlnUL - refST->offsetUL + 1;
     // The + 1 is to account for index 0 of endAlnUL

   // Set up counters for the query and reference base
   // index
   unsigned long qryNtUL = queryST->offsetUL - 1;
   unsigned long refNtUL = refST->offsetUL- 1;

   /******************************************************\
   * Fun-01 Sec-01 Sub-02:
   *  - Variables holding the scores (only two rows)
   \******************************************************/

   long snpScoreL = 0;    // Score of a single base pair
   long scoreTopL = 0;    // Score of the top cell (ins)
   long scoreDiagnolL = 0;//holds diagnol score (match/snp)
   long scoreLeftL = 0;   // Score of the left cell (del)

   // Marks when to reset score buffer (every second row)
   char swapBuffBl = 1;
   long *scoreMatrixL = 0; // matrix to use in alignment
   long *scoreOnLP = 0;  // Score I am working on
   long *topScoreLP = 0; //Pointer to cell with last base

   /******************************************************\
   * Fun-01 Sec-01 Sub-03:
   *  - Directinol matrix variables
   \******************************************************/

   // Variables dealing with the direction matrix for
   // scoring. These are as two bit array variables
   struct twoBitAry *dirMtrx = 0; // direction on
   struct twoBitAry topDir;    // insertion direction
   struct twoBitAry leftDir;   // Deletion direction

   struct scoresStruct *qryBasesST = 0;
   struct scoresStruct *refBasesST = 0;

   struct alnMatrixStruct *matrixST =
     calloc(1, sizeof(struct alnMatrixStruct));

   // Matrix scan variables

   char fileNameCStr[1024];
   char *tmpCStr = 0;
   FILE *outFILE = 0; // For matrix scan

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(matrixST == 0) return 0;

   initAlnMatrixST(matrixST);

   if(setST->matrixScanBl & 1)
   { // If I am doing a matrix scan
     strncpy(
       fileNameCStr,
       prefixCStr,
       setST->lenFileNameUS - 1
     );

     tmpCStr = fileNameCStr;

     while(*tmpCStr != '\0') ++tmpCStr;

     strncpy(
       tmpCStr,
       "--matrix-scan.aln",
       setST->lenFileNameUS - strlen(fileNameCStr) - 1
     );

     outFILE = fopen(fileNameCStr, "w");

     if(outFILE == 0)
     { // If I was unable to open the output file
       free(matrixST);
       return 0;
     } // If I was unable to open the output file
   } // If I am doing a matrix scan

   scoreMatrixL =
       calloc(2 * (lenRefUL + 1), sizeof(unsigned long));
       // I need two rows to keep track of the scores (2x)
       // lenRefUL + 1 is to account for insertion zero
       // query column
   if(scoreMatrixL == 0)
   { // If I had a memory error
     if(outFILE != 0) fclose(outFILE);
     free(matrixST);
     return 0;
   } // If I had a memory error

   dirMtrx =
     makeTwoBitArry((lenRefUL + 1) * (lenQueryUL + 1), 0);
     // Calls calloc and adds an extra element at end
       // lenRefUL + 1 accounts for insertion reference row
       // lenQeurI + 1 accounts for insertion query column

   if(dirMtrx == 0)
   { // If I do not have a direction matrix for each cell
     if(outFILE != 0) fclose(outFILE);
     free(matrixST);
     free(scoreMatrixL);
     return 0;
   } // If I do not have a direction matrix for each cell

   matrixST->dirMatrixST = dirMtrx;

   switch(setST->refQueryScanBl)
   { // Switch: Check if keeping many best scores
     case 0: break;

     case 1:
       // Make struct array for every base in the reference
       matrixST->refBasesST =
         calloc(lenRefUL, sizeof(struct scoresStruct));

       if(matrixST->refBasesST == 0)
       { // If had memory error
         if(outFILE != 0) fclose(outFILE);
         freeAlnMatrixST(matrixST);
         return 0;
       } // If had memory error

       // Make struct array for every base in the query
       matrixST->qryBasesST =
         calloc(lenQueryUL, sizeof(struct scoresStruct));

       if(matrixST->qryBasesST == 0)
       { // If had memory error
         if(outFILE != 0) fclose(outFILE);
         freeAlnMatrixST(matrixST);
         return 0;
       } // If had memory error

       qryBasesST = matrixST->qryBasesST;
       refBasesST = matrixST->refBasesST;

       matrixST->lenRefScoresUL = lenRefUL;
       matrixST->lenQueryScoresUL = lenQueryUL;

       break;
   } // Switch: Check if keeping many best scores

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Get the first indel position
   cpTwoBitPos(dirMtrx, &topDir);
   scoreOnLP = scoreMatrixL;

   // This is here just in case the user changes
   // defMoveStop from 0
   // <= is same as lenRefUL + 1 when uiCell = 0
   for(uint32_t uiCell = 0; uiCell <= lenRefUL; ++uiCell)
   { // loop; till have initalized the first row
     changeTwoBitElm(dirMtrx, defMoveStop);

     // Move to the next cell (ref base)
     ++scoreOnLP; // Already set to 0 by calloc
     twoBitAryMoveToNextElm(dirMtrx);
   } // loop; till have initalized the first row

   // Get the left (deletion0 direction positioned
   cpTwoBitPos(dirMtrx, &leftDir);
   twoBitAryMoveBackOneElm(&leftDir);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Fill the matrix with scores
   ^  o fun-01 sec-04 sub-01:
   ^    - Fill in the indel column
   ^  o fun-01 sec-04 sub-02:
   ^    - Get scores for ins,dels,matchs
   ^  o fun-01 sec-04 sub-03:
   ^    - Determine if have a best score
   ^  o fun-01 sec-4 sub-06:
   ^    - Move to the next reference base
   ^  o fun-01 sec-04 sub-04:
   ^     - Check if storing scores for multiple bases
   ^  o fun-01 sec-04 sub-05:
   ^     - Do matrix scan on previous diagnol
   ^  o fun-01 sec-04 sub-07:
   ^     - multi aligment printout, check last cell in row
   ^  o fun-01 sec-04 sub-08:
   ^     - Matrix scan, check last base in row
   ^  o fun-01 sec-04 sub-09:
   ^     - Prepare for the next row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Fill in the indel column
   \******************************************************/

   // Marks last cell in score matrix (just to make sure)
   topScoreLP = scoreMatrixL; // Set up the last base
   swapBuffBl = 1; // Swap buffers when finsh every 2nd row

   tmpQueryCStr = queryStartCStr;
   tmpRefCStr = refStartCStr;

   qryNtUL = queryST->offsetUL;
     // Not in switch, because not in loop

   // Starting on the first sequence row
   while(*tmpQueryCStr != '\0')
   { // loop; compare query base against all ref bases
     *scoreOnLP = 0;
     changeTwoBitElm(dirMtrx, defMoveStop);

     // Move to the first base comparison
     ++scoreOnLP;
     ++topScoreLP;

     // Get the directions positioned
     twoBitAryMoveToNextElm(dirMtrx);
     twoBitAryMoveToNextElm(&topDir);
     twoBitAryMoveToNextElm(&leftDir);

     tmpRefCStr = refStartCStr;

     switch(setST->multiBaseWaterBl)
     { // switch: check if keeping best score for each base
       case 1: refNtUL = refST->offsetUL;
     } // switch: check if keeping best score for each base

     /**************************************************\
     * Fun-01 Sec-04 Sub-02:
     *  - Get scores for ins, del, match
     \**************************************************/

     // First reference bases column
     while(*tmpRefCStr != '\0')
     { // loop; compare one query to one reference base
       snpScoreL =
          getBasePairScore(
              tmpQueryCStr,
              tmpRefCStr,
              setST
       ); // Find the score for the two base pairs

       // Find the score for diagnol cell (snp/match)
       scoreDiagnolL = *(topScoreLP - 1) + snpScoreL;

       scoreTopL =
           getIndelScore(
               &topDir, // Direction of last cell
               setST,
               topScoreLP
       ); // Get the score for an insertion

       scoreLeftL =
         getIndelScore(
           &leftDir,      // direction of previous cell
           setST,        // gap penalties
           scoreOnLP - 1  //Score of last base
       ); // Get the score for an insertion

       updateDirScoreWaterSingle(
         dirMtrx,
         setST,    // preference for scoring
         &scoreTopL,     // Score for insertion
         &scoreDiagnolL, // Score for match/snp
         &scoreLeftL,    // score for deletion
         scoreOnLP
       ); // Update the scores

     /**************************************************\
     * Fun-01 Sec-04 Sub-03:
     *  - Determine if have a best score
     \**************************************************/

       if(*scoreOnLP >=
          matrixST->bestScoreST.scoreL
       ) { // Else if have a new best score
           matrixST->bestScoreST.scoreL =
             *scoreOnLP;

           // Get the index of the best score

           matrixST->bestScoreST.indexUL =
             getTwoBitAryIndex(dirMtrx);
       } // Else if have a new best score

       /**************************************************\
       * Fun-01 Sec-04 Sub-04:
       *  - Check if storing scores for multiple bases
       \**************************************************/

       else
         switch(setST->multiBaseWaterBl)
         { // Switch: Check if keeping many best scores
           case 0: break;

           case 1:
           // Case 1: Keeping best score for each base
             /* What a 2x2 martix would look like
               +-++   +-++
               +b*|<--+b"| b* is the base pair to check
               +--+   +--+   if is extending.
                ^ ^     ^  b" is the base pair to check
                : :     :    after completing the row.
                : +---+ :    The idea is all b"'s turn into
                :     : :    b*'s, except for the last b".
               ++-+   +-++
               |bp|<--+bp|
               +--+   +--+
             */

             if(*(topScoreLP - 1) > setST->minScoreUI)
             { // If the previous base was worth checking
               switch(setST->refQueryScanBl)
               { // Switch; check if doing ref/query scan
                 case 0: break;

                 case 1:
                 // Case: doing a reference/qury scan
                   if(
                    getTwoBitAryElm(&leftDir) == defMoveUp
                   ){ // If the left cell continues path

                     if(
                       *(scoreOnLP-1) >= *(topScoreLP-1) &&
                       *(scoreOnLP-1) >= setST->minScoreUI
                       ) break; // If path has higher score
                   } // If the left cell continues path

                   if(
                     getTwoBitAryElm(&topDir) ==defMoveLeft
                   ){// If top cell continues the path
                     if(
                       *topScoreLP >= setST->minScoreUI &&
                       *(topScoreLP) >= *(topScoreLP - 1)
                     ) break; // If Path has a high score
                   }// If top cell continues the path

                   if(
                     getTwoBitAryElm(dirMtrx) ==
                     defMoveDiagnol
                   ){ // If the diagnol cell continues path
                     if(
                       *scoreOnLP >= setST->minScoreUI &&
                       *scoreOnLP >= *(topScoreLP - 1)
                     ) break; // If path has higher score
                   } // If the diagnol cell continues path

                   if(
                     *scoreOnLP
                       >= (qryBasesST + qryNtUL -1)->scoreL
                   ){ // If new best score for last query
                     (qryBasesST + qryNtUL - 1)->scoreL =
                       *(topScoreLP - 1);

                     (qryBasesST + qryNtUL - 1)->indexUL =
                       getTwoBitAryIndex(&leftDir);
                   } // If new best score for last query

                   else if(
                     *scoreOnLP
                       >= (refBasesST + refNtUL)->scoreL
                   ){ // else If new best score for ref
                      (refBasesST + refNtUL)->scoreL =
                       *(topScoreLP - 1);

                     (refBasesST + refNtUL)->indexUL =
                       getTwoBitAryIndex(&leftDir);
                   } // else If new best score for ref

                   break;
                 // Case: doing a reference/qury scan
               } // Switch; check if doing ref/query scan

                 /****************************************\
                 * Fun-01 Sec-04 Sub-05:
                 *  - Do matrix scan on previous diagnol
                 \****************************************/

               switch(setST->matrixScanBl)
               { // Switch; check if doing matrix scan
                 case 0: break;
                 case 1:
                 // Case: doing a matrix scan
                   // Check if This path is finished
                   if(
                    getTwoBitAryElm(&leftDir) == defMoveUp
                   ) if(*(scoreOnLP-1)>=setST->minScoreUI
                     ) break; // Path no finished

                   if(
                     getTwoBitAryElm(&topDir) ==defMoveLeft
                   ) if(*topScoreLP >=setST->minScoreUI)
                       break; // Path not finished yet

                   if(
                     getTwoBitAryElm(dirMtrx) ==
                     defMoveDiagnol
                   ) if(*scoreOnLP >= setST->minScoreUI) 
                       break; // path not finished yet

                   twoBitAryMoveBackOneElm(&topDir);

                   printMatrixCig(
                     outFILE,
                     &topDir,
                     lenRefUL,  // Length of reference
                     *(topScoreLP - 1)
                   );

                   twoBitAryMoveToNextElm(&topDir);

                   break;
                 // Case: doing a matrix scan
               } // Switch; check if doing matrix scan
             } // If the previous base was worth checking

             ++refNtUL;
             break;
           // Case 1: Keeping best score for each base
         } // Switch: Check if keeping many best scores

       /***********************************************\
       * Fun-01 Sec-04 Sub-06:
       *  - Move to next reference base
       \***********************************************/
     
       // Move to the next cell to score
       ++scoreOnLP; // next comparison for query base
       ++topScoreLP; // Move to next element
       ++tmpRefCStr;   // Move to next reference base

       // Move to the next base pair to score
       twoBitAryMoveToNextElm(dirMtrx);
       twoBitAryMoveToNextElm(&topDir);
       twoBitAryMoveToNextElm(&leftDir);
     } // loop; compare one query to one reference base

     /****************************************************\
     * Fun-01 Sec-04 Sub-07:
     *  - multi aligment printout, check last cell in row
     \****************************************************/

     switch(setST->multiBaseWaterBl)
     { // Switch: Check if keeping many best scores
       case 0: break;

       case 1:
       // Case 1: Ceck if keeping last edge cell

         if(*(topScoreLP - 1) > setST->minScoreUI)
         { // If the previous base was worth checking
          switch(setST->refQueryScanBl)
          { // Switch; check if doing ref/query scan
             case 0: break;

             case 1:
             // Case: doing a ref/query scan
               if(getTwoBitAryElm(dirMtrx) == defMoveUp)
               { // If the diagnol cell continues path
                 if(
                   *scoreOnLP >= setST->minScoreUI &&
                   *scoreOnLP >= *topScoreLP
                 ) break; // If path has higher score
               } // If the diagnol cell continues path

               if(
                 *scoreOnLP
                   >= (qryBasesST + qryNtUL - 1)->scoreL
               ){ // If new best score for last query
                 (qryBasesST + qryNtUL - 1)->scoreL =
                   *(topScoreLP);

                 (qryBasesST + qryNtUL - 1)->indexUL =
                   getTwoBitAryIndex(&topDir);
               } // If new best score for last query

               else if(
                 *scoreOnLP
                   >= (refBasesST + refNtUL)->scoreL
               ){ // else If new best score for ref
                  (refBasesST + refNtUL)->scoreL =
                   *topScoreLP;

                 (refBasesST + refNtUL)->indexUL =
                   getTwoBitAryIndex(&topDir);
               } // else If new best score for ref

               break;
             // Case: doing a ref/query scan
          } // Switch; check if doing ref/query scan

             /********************************************\
             * Fun-01 Sec-04 Sub-08:
             *  - Matrix scan, check last base in row
             \********************************************/

           switch(setST->matrixScanBl)
           { // Switch; check if doing matrix scan
             case 0: break;
             case 1:
             // Case: doing a matrix scan
               // Check if This path is finished
               if(
                getTwoBitAryElm(&leftDir) == defMoveUp
               ) if(*(scoreOnLP-1)>=setST->minScoreUI
                 ) break; // Path no finished

               if(
                 getTwoBitAryElm(&topDir) ==defMoveLeft
               ) if(*topScoreLP >=setST->minScoreUI)
                   break; // Path not finished yet

               if(
                 getTwoBitAryElm(dirMtrx) ==
                 defMoveDiagnol
               ) if(*scoreOnLP >= setST->minScoreUI) 
                   break; // path not finished yet

               twoBitAryMoveBackOneElm(&topDir);

               printMatrixCig(
                 outFILE,
                 &topDir,
                 lenRefUL,  // Length of reference
                 *topScoreLP
               );

               twoBitAryMoveToNextElm(&topDir);

               break;
             // Case: doing a matrix scan
           } // Switch; check if doing matrix scan
         } // If the previous base was worth checking

         ++qryNtUL;
         break;
       // Case 1: Keeping best score for each base
     } // Switch: Check if keeping many best scores

     /****************************************************\
     * Fun-01 Sec-04 Sub-09:
     *  - Prepare for the next row
     \****************************************************/

     if(swapBuffBl & 1)
     { // If need to swap the buffers
         scoreOnLP = scoreMatrixL; // Restart scoring
         swapBuffBl = 0;
     } // If need to swap the buffers

     else
     { // Else need to reset the score part
         topScoreLP = scoreMatrixL; // on first row
         swapBuffBl = 1;
     } // Else need to reset the score part

     ++tmpQueryCStr; // Move to the next query base
   } // loop; compare query base against all ref bases

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   ^  o fun-01 sec-05 sub-01:
   ^    - Handle final row of multibase output
   ^  o fun-01 sec-05 sub-02:
   ^    - Finsh last row for multibase output
   ^  o fun-01 sec-05 sub-03:
   ^    - Finish the last row for the matrix scan
   ^  o fun-01 sec-05 sub-04:
   ^    - Final clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-05 Sub-01:
   *  - Handle final row of multibase output
   \******************************************************/

   switch(setST->multiBaseWaterBl)
   { // Switch: Check if keeping many best scores
     case 0: break;

     case 1:
     // Case 1: Ceck if keeping last edge cell

       refNtUL = refST->offsetUL;
       qryNtUL = qryNtUL;
       tmpRefCStr = refStartCStr;

       while(*tmpRefCStr != '\0')
       { // While I have scores to check
         if(getTwoBitAryIndex(&topDir) == defMoveLeft)
           continue; // deletions just iead to first base

         if(*topScoreLP > setST->minScoreUI)
         { // If have a base worth printing out

           /**********************************************\
           * Fun-01 Sec-05 Sub-02:
           *  - Finsh last row for multibase output
           \**********************************************/
          switch(setST->refQueryScanBl)
          { // Switch; check if doing ref/query scan
             case 0: break;

             case 1:
             // Case: Not doing a matrix scan
               if(
                 *scoreOnLP >= (qryBasesST+qryNtUL)->scoreL
               ){ // If new best score for last query
                 (qryBasesST+qryNtUL)->scoreL =*topScoreLP;

                 (qryBasesST + qryNtUL)->indexUL =
                   getTwoBitAryIndex(&topDir);
               } // If new best score for last query

               else if(
                 *scoreOnLP >= (refBasesST+refNtUL)->scoreL
               ){ // else If new best score for ref
                  (refBasesST + refNtUL)->scoreL =
                   *topScoreLP;

                 (refBasesST + refNtUL)->indexUL =
                   getTwoBitAryIndex(&topDir);
               } // else If new best score for ref

               break;
             // Case: Not doing a matrix scan
          } // Switch; check if doing ref/query scan

             /********************************************\
             * Fun-01 Sec-05 Sub-03:
             *  - Finsh the last row for the matrix scan
             \********************************************/

           switch(setST->matrixScanBl)
           { // Switch; check if doing matrix scan
             case 0: break;

             case 1:
             // Case: doing a matrix scan
               printMatrixCig(
                 outFILE,
                 &topDir,
                 lenRefUL,  // Length of reference
                 *topScoreLP
               );

               break;
             // Case: doing a matrix scan
           } // Switch; check if doing matrix scan
         } // If have a base worth printing out

         ++tmpRefCStr;
       } // While I have scores to check
   } // Switch: Check if keeping many best scores
       
   /******************************************************\
   * Fun-01 Sec-05 Sub-04:
   *  - Final clean up
   \******************************************************/

   // Move back to the lower right conor cell
   twoBitAryMoveBackOneElm(dirMtrx); // not needed
   free(scoreMatrixL);
   if(outFILE != 0) fclose(outFILE);

   return matrixST;
} // WatermanAln

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints the path of the input index in dirST
\*-------------------------------------------------------*/
void printMatrixCig(
  FILE *outFILE,            // File to print to
  struct twoBitAry *dirST,  // Has index to print
  unsigned long lenRefUL,   // Length of the reference
  long scoreL               // Score of the path
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printMatrixCig
   '  - Prints out a cigar for an single path in a
   '    direction matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   uint8_t lastDirUC = defMoveStop;
   char lastDirCigC = 0;
   uint32_t numDirUI = 0;
   uint8_t curDirUC = 0;

   struct twoBitAry matrixST;

   // As index 1
   unsigned long qryEndUL =
     (getTwoBitAryIndex(dirST) / (lenRefUL + 1));

   unsigned long refEndUL =
     (getTwoBitAryIndex(dirST) % (lenRefUL + 1));

   cpTwoBitPos(dirST, &matrixST);
   curDirUC = getTwoBitAryElm(&matrixST);

   // Print out the score, position of ending query base,
   // & position of ending reference base
   fprintf(
     outFILE,
     "%ld\t%lu\t%lu\t",
     scoreL,
     qryEndUL,
     refEndUL
   );
   
   goto initializePrint;
   
   while(curDirUC != defMoveStop)
   { // While I have a cigar to build
     if(curDirUC != lastDirUC && lastDirUC != defMoveStop)
     { // If I need to print out the last direction
       if(numDirUI > 1)
         fprintf(outFILE, "%u%c", numDirUI, lastDirCigC);

       else fprintf(outFILE, "%c", lastDirCigC);

       initializePrint:

       lastDirUC = defMoveStop;
       numDirUI = 0;
       lastDirUC = curDirUC;

       switch(curDirUC)
       { // Switch: find the cigar symbol
         case defMoveUp:
           lastDirCigC = 'I';
           break;

         case defMoveDiagnol:
           lastDirCigC = 'X';
           break;

         case defMoveLeft:
           lastDirCigC = 'D';
           break;
       } // Switch: find the cigar symbol
     } // If I need to print out the last direction

     switch(curDirUC)
     { // Switch: Check which way to move
       case defMoveUp:
         twoBitAryMoveBackXElm(&matrixST, lenRefUL + 1);
         break;

       case defMoveDiagnol:
         twoBitAryMoveBackXElm(&matrixST, lenRefUL + 2);
         break;

       case defMoveLeft:
         twoBitAryMoveBackOneElm(&matrixST);
         break;
     } // Switch: Check which way to move

     ++numDirUI;
     curDirUC = getTwoBitAryElm(&matrixST);
   } // While I have a cigar to build

   // Print out the last score
   if(numDirUI > 1)
     fprintf(outFILE, "%u%c", numDirUI, lastDirCigC);

   else fprintf(outFILE, "%c", lastDirCigC);

   qryEndUL = (getTwoBitAryIndex(&matrixST) /(lenRefUL+1));
   refEndUL = (getTwoBitAryIndex(&matrixST) %(lenRefUL+1));

   // Print the starting position of query and reference
   fprintf(outFILE, "\t%lu\t%lu\n", qryEndUL, refEndUL);

   return;
} // printMatrixCig

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints out all saved alternative aliginments as
|      cigars to prefix--alt.aln
|  - Returns
|    o 1 for success
|    o 2 for file error
|    o 64 for memory error
\--------------------------------------------------------*/
unsigned char printAltWaterAlns(
  struct alnMatrixStruct *alnMtxST, // Has matrix and scores
  struct seqStruct *queryST, //query: sequence & seq length
  struct seqStruct *refST,   // ref: sequence & seq length
  struct alnSet *setST,   // Settings for the alignment
  char *prefxCStr            // Prefix of file to write to
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: printAltWaterAlns
   '  - Prints out all saved alternatives alignments
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Open the output file
   '  o fun-04 sec-03:
   '    o Print out the reference alignments
   '  o fun-04 sec-04:
   '    - Print out the query aligments
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-01:
  ^  - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  unsigned short lenFileNameUS = 1024;
  char fileNameCStr[lenFileNameUS];
  char *tmpCStr = 0;

  unsigned long lenRefUL = refST->lenSeqUL-refST->offsetUL;

  struct scoresStruct *scoreST = 0;
  struct twoBitAry dirOn;
   
  FILE *outFILE = 0;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-02:
  ^  - Open the output file
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  strncpy(fileNameCStr, prefxCStr, lenFileNameUS - 1);
  tmpCStr = fileNameCStr;

  while(*tmpCStr != '\0') ++tmpCStr;

  strncpy(
    tmpCStr,
    "--alt.aln",
    lenFileNameUS - strlen(fileNameCStr) - 1
  ); // Add the ending to the file name

  outFILE = fopen(fileNameCStr, "w");

  if(outFILE == 0) return 2;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-03:
  ^  - Print out alternative alignments as cigars
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  cpTwoBitPos(alnMtxST->dirMatrixST, &dirOn);
  scoreST = alnMtxST->refBasesST;

  for(
    unsigned long ulRefBase = 0;
    ulRefBase < alnMtxST->lenRefScoresUL;
    ++ulRefBase
  ){ // For all reference bases in the alignment

    // Check if even need to print out any reference alns
    if(scoreST->scoreL < setST->minScoreUI)
      goto nextRefScore;

    moveXElmFromStart(&dirOn, scoreST->indexUL);

    printMatrixCig(
      outFILE,
      &dirOn,
      lenRefUL,
      scoreST->scoreL
    ); // Prit out the cigar entry

     nextRefScore:
    ++scoreST;  // Move to the next entry
  } // For all reference bases in the alignment
     
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Print out the query alignments
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  scoreST = alnMtxST->qryBasesST;

  for(
    unsigned long ulQryBase = 0;
    ulQryBase < alnMtxST->lenQueryScoresUL;
    ++ulQryBase
  ){ // For all reference bases in the alignment

    // Check if even need to print out any reference alns
    if(scoreST->scoreL < setST->minScoreUI)
      goto nextQueryScore;

    moveXElmFromStart(&dirOn, scoreST->indexUL);

    printMatrixCig(
      outFILE,
      &dirOn,
      lenRefUL,
      scoreST->scoreL
    ); // Prit out the cigar entry

    nextQueryScore:
    ++scoreST;  // Move to the next entry
  } // For all reference bases in the alignment

  fclose(outFILE);

  return 1;
} // printAltWaterAlns

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o scoreOnL to hold the best score
|    o dirOnUC to hold the best direction
\--------------------------------------------------------*/
void updateDirScoreWaterSingle(
    struct twoBitAry *dirOnST, // has matrix cell to update
    struct alnSet *alnSetST,
      // Has preference for selecting equal scores
    long *scoreTopL,     // Score for an insertion
    long *scoreDiagnolL, // Score for an match/snp
    long *scoreLeftL,    // The score for an deletion
    long *scoreOnL       // Score to update
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: updateDirScoreWaterSingle
   '  - Picks the best score and direction for the current
   '    base pairs being compared in a Waterman Smith
   '    alignment
   '  o fun-05 sec-1: Matches->insertions->deletions
   '  o fun-05 sec-2: Matches->deletions->insertions
   '  o fun-05 sec-3: Insertions->matches->deletions
   '  o fun-05 sec-4: Deletions->matches->insertions
   '  o fun-05 sec-5: Insertions->deletions->matches
   '  o fun-05 sec-6: Deletions->insertions->matches
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-1: Matches and then insertions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   switch(alnSetST->diagnolPriorityC)
   { // Switch; get an snp/match priority
       case 0:  // Top priority
       // Case; bases or matches are top priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 1:
               // Case: priority matches/snps then ins
                   if(*scoreDiagnolL >= *scoreTopL)
                   { // If diagnol beats insertion
                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );

                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveDiagnol
                           );

                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats deletions

                       else
                       { // Else deletion is the best score
                           if(*scoreLeftL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );

                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveLeft
                           );
                           *scoreOnL = *scoreLeftL;
                           return;
                       } // Else the deletion is the best
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else the insertion is best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(dirOnST, defMoveUp);
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else the insertion is best score

                   else
                   { // Else the deletion is best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveLeft
                      );
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else the deletion is best score
                       
                   return;
               // Case: priority matches/snps then ins

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-05 Sec-2: Matches->dels->ins
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority matches/snps then dels
                   if(*scoreDiagnolL >= *scoreLeftL)
                   { // If diagnol beats deletion
                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );

                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveDiagnol
                           );
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats insertions

                       else
                       { // Else insertion is best score
                          if(*scoreTopL <= 0)
                          { // If need to make a stop
                              changeTwoBitElm(
                                dirOnST,
                                defMoveStop
                              );
                              *scoreOnL = 0;
                              return;
                          } // If need to make a stop

                          changeTwoBitElm(
                             dirOnST,
                             defMoveUp
                           );
                          *scoreOnL = *scoreTopL;
                          return;
                       } // Else insertion is best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else the deletion is best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                            dirOnST,
                            defMoveStop
                          );

                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveLeft
                      );
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else deletion is best score

                   else
                   { // Else insertion is best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveUp
                      );
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else insertion is best score
                       
                   return;
               // Case: priority matches/snps then dels
           } // Priority for insertions
       // Case; bases or matches are top priority

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-05 Sec-3: Insertions->matches->deletions
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 1:
       // Case; bases or matches are second priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority ins then matches/snps
                   if(*scoreDiagnolL > *scoreTopL)
                   { // If diagnol beats insertion
                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveDiagnol
                           );
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats deletions

                       else
                       { // Else deletion is best score
                           if(*scoreLeftL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveLeft
                           );
                           *scoreOnL = *scoreLeftL;
                           return;
                       } // Else deletion is best score
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else insertion is best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveUp
                      );
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else insertion is best score

                   else
                   { // Else deletion is best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveLeft
                      );
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else deletion is best score
                       
                   return;
               // Case: priority matches/snps then ins

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-05 Sec-4: Dels->matches->ins
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority dels then matches/snps
                   if(*scoreDiagnolL > *scoreLeftL)
                   { // If diagnol beats deletion

                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveDiagnol
                           );
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats insertions

                       else
                       { // Else ins is the best score
                           if(*scoreTopL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                  dirOnST,
                                  defMoveStop
                                );
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveUp
                           );
                           *scoreOnL = *scoreTopL;
                           return;
                       } // Else insertion is best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else del is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveLeft
                      );
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else deletion is best score

                   else
                   { // Else insertion is best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveUp
                      );
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else insertion is best score
                       
                   return;
               // Case: priority matches/snps then del
           } // Priority for insertions

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-05 Sec-5: Insertions->deletions->matches
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 2:
       // Case; bases or matches are last priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority insertions then deletions
                   if(*scoreTopL >= *scoreLeftL)
                   { // If diagnol beats insertion
                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol is highest score
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveDiagnol
                           );
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol is highest score

                       else
                       { // Else insertion is best score
                           if(*scoreTopL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );

                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveUp
                           );
                           *scoreOnL = *scoreTopL;
                           return;
                       } // Else deletion is best score
                   } // If diagnol beats insertion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveLeft
                      );
                      *scoreOnL = *scoreLeftL;
                   } // Else deletion is the best score

                   else
                   { // Else match/snp is the best score
                      if(*scoreDiagnolL <= 0)
                      { // If need to make a stopping point
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      changeTwoBitElm(
                        dirOnST,
                        defMoveDiagnol
                      );
                      *scoreOnL = *scoreDiagnolL;
                      return;
                   } // Else match/snp is the best score
                       
                   return;
               // Case: priority insertions then deletions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-05 Sec-6: Dels->insertions->matches
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority deletions then insertions
                   if(*scoreTopL > *scoreLeftL)
                   { // If insertion beats deletion
                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol beats insertions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                  dirOnST,
                                  defMoveStop
                                );
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveDiagnol
                           );
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats insertions

                       else
                       { // Else insertion is best score
                           if(*scoreTopL <= 0)
                           { // If need to make a stop
                               changeTwoBitElm(
                                 dirOnST,
                                 defMoveStop
                               );
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stop

                           changeTwoBitElm(
                             dirOnST,
                             defMoveUp
                           );
                           *scoreOnL = *scoreTopL;
                           return;
                       } // Else insertion is best score
                   } // If insertion beats deletion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(dirOnST,defMoveLeft);
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else deletion is the best score

                   else
                   { // Else match/snp is the best score
                      if(*scoreDiagnolL <= 0)
                      { // If need to make a stop
                          changeTwoBitElm(
                             dirOnST,
                             defMoveStop
                           );
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stop

                      changeTwoBitElm(
                        dirOnST,
                        defMoveDiagnol
                      );
                      *scoreOnL = *scoreDiagnolL;
                      return;
                   } // Else match/snp is the best score
                       
                   return;
               // Case: priority insertions then deletions
           } // Priority for insertions
       // Case; bases or matches are last priority
   } // Switch; get an snp/match priority

   return;
} // updateDirScoreWaterSingle
