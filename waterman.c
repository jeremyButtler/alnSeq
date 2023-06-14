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
#########################################################*/

#include "waterman.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 WatermanSmithAln:
'     o Perform a Waterman Smith alignment on input
'       sequences
'  o fun-02 addBestBaseScore:
'    - Adds a score and index to the kept scores list
'  o fun-03 printAltWaterAlns:
'    - Prints out the best aligment and the saved
'       alterantive alignments  (best alignment for each
'       base) to a file
'  - fun-04 updateDirScoreWaterSingle:
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
    struct seqStruct *queryST,
      // query sequence, length, & bounds for alignment
    struct seqStruct *refST,  
      // reference sequence, length, & bounds for alignment
    struct alnSet *settings // Settings for the alignment
    // *startI and *endI paramaters should be index 1
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
   char *refStartCStr = refST->seqCStr + refST->offsetUI;

   char *queryStartCStr =
        queryST->seqCStr + queryST->offsetUI;

   char *tmpQueryCStr = queryStartCStr;
   char *tmpRefCStr = refStartCStr;

   // Find the length of the reference and query
   unsigned long lenQueryUL =
       queryST->endAlnUI - queryST->offsetUI;

   unsigned long lenRefUL =
       refST->endAlnUI - refST->offsetUI;

   // Set up counters for the query and reference base
   // index
   unsigned long queryNtOnUL = queryST->offsetUI - 1;
   unsigned long refNtOnUL = refST->offsetUI- 1;

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
   long *scoreOnLPtr = 0;  // Score I am working on
   long *lastBaseLPtr = 0; //Pointer to cell with last base

   /******************************************************\
   * Fun-01 Sec-01 Sub-03:
   *  - Directinol matrix variables
   \******************************************************/

   // Variables dealing with the direction matrix for
   // scoring. These are as two bit array variables
   struct twoBitAry *dirMatrix = 0; // direction on
   struct twoBitAry topDir;    // insertion direction
   struct twoBitAry leftDir;   // Deletion direction

   struct scoresStruct **lastQueryScoreST = 0;
   struct scoresStruct **lastRefScoreST = 0;

   struct alnMatrixStruct *retMtxST =
     calloc(1, sizeof(struct alnMatrixStruct));

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(retMtxST == 0) return 0;

   initAlnMatrixST(retMtxST);

   scoreMatrixL =
       calloc(2 * (lenRefUL + 1), sizeof(unsigned long));
       // I need two rows to keep track of the scores (2x)
       // lenRefUL + 1 is to account for insertion zero
       // query column
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

   if(dirMatrix == 0)
   { // If I do not have a direction matrix for each cell
     free(retMtxST);
     free(scoreMatrixL);
     return 0;
   } // If I do not have a direction matrix for each cell

   retMtxST->dirMatrixST = dirMatrix;

   switch(settings->multiBaseWaterBl)
   { // Switch: Check if keeping many best scores
     case 0: break;

     case 1:
       // Make struct array for every base in the reference
       retMtxST->refScoresST =
         calloc(lenRefUL, sizeof(struct scoresStruct));

       if(retMtxST->refScoresST == 0)
       { // If had memory error
         freeAlnMatrixST(retMtxST);
         return 0;
       } // If had memory error

       // Make struct array for every base in the query
       retMtxST->queryScoresST =
         calloc(lenQueryUL, sizeof(struct scoresStruct));

       if(retMtxST->queryScoresST == 0)
       { // If had memory error
         freeAlnMatrixST(retMtxST);
         return 0;
       } // If had memory error

       retMtxST->lenRefScoresUL = lenRefUL;
       retMtxST->lenQueryScoresUL = lenQueryUL;

       lastQueryScoreST =
         calloc(lenQueryUL, sizeof(struct scoresStruct));

       if(lastQueryScoreST == 0)
       { // If had memory error
         freeAlnMatrixST(retMtxST);
         return 0;
       } // If had memory error


       lastRefScoreST =
         calloc(lenRefUL, sizeof(struct scoresStruct));

       if(lastRefScoreST == 0)
       { // If had memory error
         free(lastQueryScoreST);
         freeAlnMatrixST(retMtxST);
         return 0;
       } // If had memory error

       break;
   } // Switch: Check if keeping many best scores

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Calloc already set everything to done (defMoveStop=0)
   scoreOnLPtr = scoreMatrixL + lenRefUL + 1;
       // Move to the first column in next row

   // Get the first indel position
   cpTwoBitPos(dirMatrix, &topDir);

   // Move to start of the first sequence indel column
   twoBitAryMoveForXElm(dirMatrix, lenRefUL + 1);

   // Get the left (deletion0 direction positioned
   cpTwoBitPos(dirMatrix, &leftDir);
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
   ^  o fun-01 sec-4 sub-04:
   ^    - Move to the next base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-04 Sub-01:
   *  - Fill in the indel column
   \******************************************************/

   // Marks last cell in score matrix (just to make sure)
   lastBaseLPtr = scoreMatrixL; // Set up the last base
   swapBuffBl = 1; // Swap buffers when finsh every 2nd row

   tmpQueryCStr = queryStartCStr;
   tmpRefCStr = refStartCStr;

   // Starting on the first sequence row
   while(*tmpQueryCStr != '\0')
   { // loop; compare query base against all ref bases
     *scoreOnLPtr = 0;
     changeTwoBitElm(dirMatrix, defMoveStop);

     // Move to the first base comparison
     ++scoreOnLPtr;
     ++lastBaseLPtr;

     // Get the directions positioned
     twoBitAryMoveToNextElm(dirMatrix);
     twoBitAryMoveToNextElm(&topDir);
     twoBitAryMoveToNextElm(&leftDir);

     tmpRefCStr = refStartCStr;

     switch(settings->multiBaseWaterBl)
     { // switch: check if keeping best score for each base
       case 1: refNtOnUL = refST->offsetUI - 1;
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
              settings
       ); // Find the score for the two base pairs

       // Find the score for diagnol cell (snp/match)
       scoreDiagnolL = *(lastBaseLPtr - 1) + snpScoreL;

       scoreTopL =
           getIndelScore(
               &topDir, // Direction of last cell
               settings,
               lastBaseLPtr
       ); // Get the score for an insertion

       scoreLeftL =
         getIndelScore(
           &leftDir,      // direction of previous cell
           settings,        // gap penalties
           scoreOnLPtr - 1  //Score of last base
       ); // Get the score for an insertion

       updateDirScoreWaterSingle(
         dirMatrix,
         settings,    // preference for scoring
         &scoreTopL,     // Score for insertion
         &scoreDiagnolL, // Score for match/snp
         &scoreLeftL,    // score for deletion
         scoreOnLPtr
       ); // Update the scores

     /**************************************************\
     * Fun-01 Sec-04 Sub-03:
     *  - Determine if have a best score
     \**************************************************/

       if(*scoreOnLPtr >=
          retMtxST->bestScoreST.scoreL
       ) { // Else if have a new best score
           retMtxST->bestScoreST.scoreL =
             *scoreOnLPtr;

           // Get the index of the best score

           retMtxST->bestScoreST.indexUL =
             getTwoBitAryIndex(dirMatrix);
       } // Else if have a new best score

       else
         switch(settings->multiBaseWaterBl)
         { // Switch: Check if keeping many best scores
           case 0: break;

           case 1:
           // Case 1: Keeping best score for each base
             addBestBaseScore(
               getTwoBitAryElm(dirMatrix),
               getTwoBitAryIndex(dirMatrix),
               *scoreOnLPtr,
               lenRefUL,
               *retMtxST->refScoresST,
               refNtOnUL,
               *retMtxST->queryScoresST,
               queryNtOnUL,
               *lastRefScoreST,
               *lastQueryScoreST 
             ); // Add the new score in

             ++refNtOnUL;
             break;
           // Case 1: Keeping best score for each base
         } // Switch: Check if keeping many best scores

       /***********************************************\
       * Fun-01 Sec-04 Sub-04:
       *  - Move to next bases
       \***********************************************/
     
       // Move to the next cell to score
       ++scoreOnLPtr; // next comparison for query base
       ++lastBaseLPtr; // Move to next element
       ++tmpRefCStr;   // Move to next reference base

       // Move to the next base pair to score
       twoBitAryMoveToNextElm(dirMatrix);
       twoBitAryMoveToNextElm(&topDir);
       twoBitAryMoveToNextElm(&leftDir);
     } // loop; compare one query to one reference base

     if(swapBuffBl & 1)
     { // If need to swap the buffers
         scoreOnLPtr = scoreMatrixL; // Restart scoring
         swapBuffBl = 0;
     } // If need to swap the buffers

     else
     { // Else need to reset the score part
         lastBaseLPtr = scoreMatrixL; // on first row
         swapBuffBl = 1;
     } // Else need to reset the score part

     switch(settings->multiBaseWaterBl)
     { // switch: check if keeping best score for each base
       case 1: ++queryNtOnUL;
     } // switch: check if keeping best score for each base

     ++tmpQueryCStr; // Move to the next query base
   } // loop; compare query base against all ref bases

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Move back to the lower right conor cell
   twoBitAryMoveBackOneElm(dirMatrix); // not needed
   free(scoreMatrixL);

   return retMtxST;
} // WatermanAln

/*--------------------------------------------------------\
| Output:
|  - Modifes:
|    o socresST to hold the new score and index if better
|      than the old score and index
|    o oldScoreST to hold the old score if it is part of
|      a different alignment
\--------------------------------------------------------*/
void addBestBaseScore(
  uint8_t dirUC,                   // Direction travled
  long indexUL,                 // new index to add in
  long scoreL,                 // new score to add in
  unsigned long lenRefUI,  // length of reference sequence
  struct scoresStruct *refScoreST,//all Kept ref scores
  unsigned long refBaseUL,       // current ref base on
  struct scoresStruct *queryScoreST,//all kept query scores
  unsigned long queryBaseUL,     // Currnetn query base on
  struct scoresStruct *oldRScoreST,// Has old ref scores
  struct scoresStruct *oldQScoreST // Has old query scores
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-01: addBestBaseScore
   '  - Adds a score and index to the kept scores list
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  struct scoresStruct *scoreOnST = 0;
  struct scoresStruct *oldScoreOnST = 0;
  
  if((refScoreST + refBaseUL)->scoreL < scoreL)
  { // If replacing the reference score
    scoreOnST = refScoreST + refBaseUL;
    oldScoreOnST = oldRScoreST + refBaseUL;
  } // If replacing the reference score

  else if((queryScoreST + queryBaseUL)->scoreL < scoreL)
  { // Else if I am replacing thw query base
    scoreOnST = queryScoreST + queryBaseUL;
    oldScoreOnST = oldQScoreST + queryBaseUL;
  } // Else if I am replacing thw query base

  else return; // Nothing to do

  switch(dirUC)
  { // Switch: check last direction
    // For stops I know there is no worry about pointing
    // to a previous alignment
    case defMoveStop:
    // Case: This is a stop direction
      oldScoreOnST->indexUL = scoreOnST->indexUL;
      oldScoreOnST->scoreL = scoreOnST->scoreL;

      scoreOnST->indexUL = indexUL;
      scoreOnST->scoreL = scoreL;

      return;
    // Case: This is a stop direction

    case defMoveLeft:
    // Case: Moved left (deletion)
      // Check if I am continuing an already recorded
      // path. If so, I have no need of keeping the last
      // reference bases stored best score or the query
      // bases score
      if((refScoreST + refBaseUL - 1)->indexUL ==indexUL-1)
        { // If I am contuning on the same path
           (refScoreST + refBaseUL - 1)->indexUL = 
             (oldRScoreST + refBaseUL - 1)->indexUL;

           (refScoreST + refBaseUL - 1)->scoreL = 
             (oldRScoreST + refBaseUL - 1)->scoreL;
        } // If I am contuning on the same path

        // Still on the same query base, but not same ref
        else if(
          (queryScoreST + queryBaseUL - 1)->indexUL
         == 
          indexUL - 1
        ) { // Else if contnuing the qury's path
           (queryScoreST + queryBaseUL)->indexUL = 
             (oldQScoreST + queryBaseUL)->indexUL;

           (queryScoreST + queryBaseUL)->scoreL = 
             (oldQScoreST + queryBaseUL)->scoreL;
        } // Else if contnuing the qury's path

        scoreOnST->indexUL = indexUL;
        scoreOnST->scoreL = scoreL;
 
        return;
    // Case: Moved left (deletion)

    case defMoveUp:
    // Case: Moved up (insertion)
      // Check if I am continuing an already recorded
      // path. If so, I have no need of keeping the last
      // reference bases stored best score or the query
      // bases score
      // On same reference, but not query
      if((
         refScoreST + refBaseUL)->indexUL
        ==
         indexUL - lenRefUI - 1
        ){ // If I am contuning on the same path
           (refScoreST + refBaseUL)->indexUL = 
             (oldRScoreST + refBaseUL)->indexUL;

           (refScoreST + refBaseUL)->scoreL = 
             (oldRScoreST + refBaseUL)->scoreL;
        } // If I am contuning on the same path

        else if(
          (queryScoreST + queryBaseUL - 1)->indexUL
         == 
          indexUL - lenRefUI - 1
        ) { // Else if contnuing the qury's path
           (queryScoreST + queryBaseUL - 1)->indexUL = 
             (oldQScoreST + queryBaseUL - 1)->indexUL;

           (queryScoreST + queryBaseUL - 1)->scoreL = 
             (oldQScoreST + queryBaseUL - 1)->scoreL;
        } // Else if contnuing the qury's path

        scoreOnST->indexUL = indexUL;
        scoreOnST->scoreL = scoreL;
 
        return;
    // Case: Moved up (insertion)

    case defMoveDiagnol:
    // Case: Moved diaginol (match/snp)
      // Check if I am continuing an already recorded
      // path. If so, I have no need of keeping the last
      // reference bases stored best score or the query
      // bases score
      // In this case both the query and reference change
      if((
         refScoreST + refBaseUL)->indexUL
        ==
         indexUL - lenRefUI - 2
        ){ // If I am contuning on the same path
           (refScoreST + refBaseUL - 1)->indexUL = 
             (oldRScoreST + refBaseUL - 1)->indexUL;

           (refScoreST + refBaseUL - 1)->scoreL = 
             (oldRScoreST + refBaseUL - 1)->scoreL;
        } // If I am contuning on the same path

        else if(
          (queryScoreST + queryBaseUL - 1)->indexUL
         == 
          indexUL - lenRefUI - 2
        ) { // Else if contnuing the qury's path
           (queryScoreST + queryBaseUL - 1)->indexUL = 
             (oldQScoreST + queryBaseUL - 1)->indexUL;

           (queryScoreST + queryBaseUL - 1)->scoreL = 
             (oldQScoreST + queryBaseUL - 1)->scoreL;
        } // Else if contnuing the qury's path

        scoreOnST->indexUL = indexUL;
        scoreOnST->scoreL = scoreL;
 
        return;
    // Case: Moved diaginol (match/snp)
  } // Switch: check last direction

  return;
} // addBestBaseScore

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints the best aligment and each alignment that is
|      kept to a separate file
|    o Best ailgment:
|      prefxCStr--Best--ref-start-end--query-start-end.aln
|    o Alignments from the reference sequence:
|      prefx--Reference--ref-start-end--query-start-end.aln
|    o Alignments from the Query sequence:
|      prefxStr--Query--ref-start-end--query-start-end.aln
|  - Returns
|    o 1 for success
|    o 2 for file error
|    o 64 for memory error
\--------------------------------------------------------*/
unsigned char printAltWaterAlns(
  struct alnMatrixStruct *alnMtxST, // Has matrix and scores
  struct seqStruct *queryST, //query: sequence & seq length
  struct seqStruct *refST,   // ref: sequence & seq length
  struct alnSet *settings,   // Settings for the alignment
  char *prefxCStr            // Prefix of file to write to
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printAltWaterAlns
   '  - Prints out the best aligment and the saved
   '     alterantive alignments  (best alignment for each
   '     base) to a file
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Print out the best alignment
   '  o fun-03 sec-03:
   '    o Sort scores to get highest so lowest scores
   '  o fun-03 sec-04:
   '    - Print out the reference scores
   '  o fun-03 sec-05:
   '    - Get the query alignments
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-01:
  ^  - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  unsigned short lenFileNameUS = 1024;
  char *queryAlnCStr = 0;
  char *refAlnCStr = 0;
  char fileNameCStr[lenFileNameUS];

  struct scoresStruct *scoreST = 0;
  struct alnStruct *alnST = 0;
  FILE *outFILE = 0;

  // I removed the needed for this. I left the code
  // that used this commented out for later use
  //struct scoresStruct *refScoreST = 0;
  //uint32_t refBaseUI = 0;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-02:
  ^  - Print out the best alignment
  ^  o fun-03 sec-02 sub-01:
  ^    - Find alignment and make alignment human readable
  ^  o fun-03 sec-02 sub-02:
  ^    - Make output file name and print out alingment
  ^  o fun-03 sec-02 sub-03:
  ^    - Free uneeded variables (clean up)
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*******************************************************\
  * Fun-03 Sec-02 Sub-01:
  *  - Find alignment and make alignment human readable
  \*******************************************************/

  // Get the alignment
  alnST =
    cnvtDirMatrixToAlnAry(
      refST,
      queryST,
      alnMtxST->dirMatrixST,
      &alnMtxST->bestScoreST,
      1      // Applying a soft mask
  );

  if(alnST == 0) return 2;
  
  // Align the reference and query bases to the alignment
  refAlnCStr = cnvtAlnAryToSeq(refST, 0, alnST);

  if(refAlnCStr == 0)
  { // If had a memroy error
    freeAlnST(alnST, 1); // No longer need
    return 64;
  } // If had a memroy error

  queryAlnCStr = cnvtAlnAryToSeq(queryST, 0, alnST);

  if(queryAlnCStr == 0)
  { // If had a memroy error
    free(refAlnCStr);
    freeAlnST(alnST, 1); // No longer need
    return 64;
  } // If had a memroy error

  // Conver the alingment codes to human readable
  alnAryToLetter(refST, queryST, alnST);

  /*******************************************************\
  * Fun-03 Sec-02 Sub-02:
  *  - Make output file name and print out alingment
  \*******************************************************/

  // Make the output file and open it
  snprintf(
    fileNameCStr,
    lenFileNameUS,
    "%s--Best--ref-%u-%u--query-%u-%u.aln",
    prefxCStr,
    alnST->refStartUI,
    alnST->refEndUI,
    alnST->queryStartUI,
    alnST->queryEndUI
  ); // Make the new file name

  outFILE = fopen(fileNameCStr, "w");

  if(outFILE == 0)
  { // If I could not open the file
    free(refAlnCStr);
    free(queryAlnCStr);
    freeAlnST(alnST, 1); // No longer need

    return 1;
  } // If I could not open the file

  // Print the alignment
  printAln(
    outFILE,
    queryST->idCStr,
    queryAlnCStr,
    refAlnCStr,
    refST->idCStr,
    alnMtxST->bestScoreST.scoreL,
    settings->lineWrapUS,
    alnST
  );

  /*******************************************************\
  * Fun-03 Sec-02 Sub-03:
  *  - Free uneeded variables (clean up)
  \*******************************************************/

  free(refAlnCStr);
  refAlnCStr = 0;

  free(queryAlnCStr);
  queryAlnCStr = 0;

  freeAlnST(alnST, 1); // No longer need
  alnST = 0;

  fclose(outFILE);

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-03:
  ^  - Sort scores to get highest so lowest scores
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /* This was for sorthing scores, currently not needed
  // Sort the reference scores (sorts > to least)
  sortScores(
    alnMtxST->refScoresST,
    refST->offsetUI,
    refST->endAlnUI
  );

  // Sort the scores for the query (sorts > to least)
  sortScores(
    alnMtxST->queryScoresST,
    queryST->offsetUI,
    queryST->endAlnUI
  );
  */

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Print out the reference scores
  ^  o fun-03 sec-04 sub-01:
  ^    - Make sure the alignment is worth printing out
  ^  o fun-03 sec-04 sub-02:
  ^    - Find alignment and make alignment human readable
  ^  o fun-03 sec-04 sub-03:
  ^    - Make output file name and print out alingment
  ^  o fun-03 sec-04 sub-04:
  ^    - Free uneeded variables (clean up)
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*******************************************************\
  * Fun-03 Sec-04 Sub-01:
  *  - Make sure the alignment is worth printing out
  \*******************************************************/
  
  scoreST = *alnMtxST->refScoresST;

  for(
    unsigned long ulRefBase = 0;
    ulRefBase < alnMtxST->lenRefScoresUL;
    ++ulRefBase
  ){ // For all reference bases in the alignment

    // No longer needed. I mades sure this could never
    // happen
    //// Check if this is the same index as the best score
    //if(scoreST->indexUL == alnMtxST->bestScoreST.indexUL)
    //{ // If on the best score (move to next score)
    //  ++scoreST; 
    //  continue;
    //} // If on the best score (move to next score)

    // Check if even need to print out any reference alns
    if(scoreST->scoreL < settings->minScoreUI)
      break; // No more socres to print out

    /*****************************************************\
    * Fun-03 Sec-04 Sub-02:
    *  - Find alignment and make alignment human readable
    \*****************************************************/

    alnST = 
      cnvtDirMatrixToAlnAry(
        refST,
        queryST,
        alnMtxST->dirMatrixST,
        scoreST,
        1      // Applying a soft mask
    );

    if(alnST == 0) return 2;

    if(alnST->numAlnBasesUI < settings->minBasesUI)
    { // If the alignment does not have enough bases
      freeAlnST(alnST, 1); // No longer need
      ++scoreST;
      continue;  // Move to the next alignmetn
    } // If the alignment does not have enough bases
    
    // Align the reference and query bases to the alignment
    refAlnCStr = cnvtAlnAryToSeq(refST, 0, alnST);
  
    if(refAlnCStr == 0)
    { // If had a memroy error
      freeAlnST(alnST, 1); // No longer need
      return 64;
    } // If had a memroy error
  
    queryAlnCStr = cnvtAlnAryToSeq(queryST, 0, alnST);
  
    if(queryAlnCStr == 0)
    { // If had a memroy error
      free(refAlnCStr);
      freeAlnST(alnST, 1); // No longer need
      return 64;
    } // If had a memroy error
  
    // Conver the alingment codes to human readable
    alnAryToLetter(refST, queryST, alnST);
  
    /*******************************************************\
    * Fun-03 Sec-04 Sub-03:
    *  - Make output file name and print out alingment
    \*******************************************************/
  
    // Make the output file and open it
    snprintf(
      fileNameCStr,
      lenFileNameUS,
      "%s--Reference--ref-%u-%u--query-%u-%u.aln",
      prefxCStr,
      alnST->refStartUI,
      alnST->refEndUI,
      alnST->queryStartUI,
      alnST->queryEndUI
    ); // Make the new file name
  
    outFILE = fopen(fileNameCStr, "w");
  
    if(outFILE == 0)
    { // If I could not open the file
      free(refAlnCStr);
      free(queryAlnCStr);
      freeAlnST(alnST, 1); // No longer need
  
      return 1;
    } // If I could not open the file
  
    // Print the alignment
    printAln(
      outFILE,
      queryST->idCStr,
      queryAlnCStr,
      refAlnCStr,
      refST->idCStr,
      alnMtxST->bestScoreST.scoreL,
      settings->lineWrapUS,
      alnST
    );
  
    /*******************************************************\
    * Fun-03 Sec-04 Sub-04:
    *  - Free uneeded variables (clean up)
    \*******************************************************/
  
    free(refAlnCStr);
    refAlnCStr = 0;
  
    free(queryAlnCStr);
    queryAlnCStr = 0;
  
    freeAlnST(alnST, 1); // No longer need
    alnST = 0;
  
    fclose(outFILE);
    ++scoreST;  // Move to the next entry
  } // For all reference bases in the alignment
     
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-05:
  ^  - Get the query alignments
  ^  o Fun-03 sec-05 sub-01:
  ^    - Set up for query aligment printing
  ^  o fun-03 sec-05 sub-02:
  ^    - Check if I have already printed this alignment
  ^  o fun-03 sec-05 sub-03:
  ^    - Find alignment and make alignment human readable
  ^  o fun-03 sec-05 sub-04:
  ^    - Make output file name and print out alingment
  ^  o fun-03 sec-05 sub-05:
  ^    - Free uneeded variables (clean up)
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*******************************************************\
  * Fun-03 Sec-05 Sub-01:
  *  - Set up for query aligment printing
  \*******************************************************/

  scoreST = *alnMtxST->queryScoresST;
  //refScoreST = *alnMtxST->refScoresST;
  //refBaseUI = refST->offsetUI;

  // This was for when I had sorted scores
  //if(scoreST->scoreL < settings->minScoreUI)
  //  return;

  for(
    uint32_t ulQueryBase = 0;
    ulQueryBase < alnMtxST->lenQueryScoresUL;
    ++ulQueryBase
  ){ // For all reference bases in the alignment

    /*****************************************************\
    * Fun-03 Sec-05 Sub-02:
    *  - Check if I have already printed this alignment
    \*****************************************************/

    /* I made sure these cases could never happen
    // Check if this is the same index as the best score
    if(scoresST->indexUL == alnMtxST->bestScore->indexUL)
    { // If on the best score (move to next score)
      ++scoreST; 
      continue;
    } // If on the best score (move to next score)

    // Check if even need to print out any reference alns
    if(scoreST->scoreL < settings->minScoreUI)
      break; // No more socres to print out

    if(refBaseUI >= refST->endAlnUI) goto afterRefCheck;

    while(refScoreST->scoreL > scoreST->scoreL)
    { // Adavnce to the next reference score
      ++refScoreST;
      ++refBaseUI;
      if(refBaseUI >= refST->endAlnUI) goto afterRefCheck;
    } // Adavnce to the next reference score

    if(refScoreST->indexUL == scoreST->indexUL)
    { // If I have already printed this alignment as ref
      ++scoreST;
      continue;
    } // If I have already printed this alignment as ref

    afterRefCheck: // If no more reference scores to check
    */

    /*****************************************************\
    * Fun-03 Sec-05 Sub-03:
    *  - Find alignment and make alignment human readable
    \*****************************************************/

    alnST = 
      cnvtDirMatrixToAlnAry(
        refST,
        queryST,
        alnMtxST->dirMatrixST,
        scoreST,
        1      // Applying a soft mask
    );

    if(alnST == 0) return 2;

    if(alnST->numAlnBasesUI < settings->minBasesUI)
    { // If the alignment does not have enough bases
      freeAlnST(alnST, 1); // No longer need
      ++scoreST;
      continue;  // Move to the next alignmetn
    } // If the alignment does not have enough bases
    
    // Align the reference and query bases to the alignment
    refAlnCStr = cnvtAlnAryToSeq(refST, 0, alnST);
  
    if(refAlnCStr == 0)
    { // If had a memroy error
      freeAlnST(alnST, 1); // No longer need
      return 64;
    } // If had a memroy error
  
    queryAlnCStr = cnvtAlnAryToSeq(queryST, 0, alnST);
  
    if(queryAlnCStr == 0)
    { // If had a memroy error
      free(refAlnCStr);
      freeAlnST(alnST, 1); // No longer need
      return 64;
    } // If had a memroy error
  
    // Conver the alingment codes to human readable
    alnAryToLetter(refST, queryST, alnST);
  
    /*****************************************************\
    * Fun-03 Sec-05 Sub-04:
    *  - Make output file name and print out alingment
    \*****************************************************/
  
    // Make the output file and open it
    snprintf(
      fileNameCStr,
      lenFileNameUS,
      "%s--Query--ref-%u-%u--query-%u-%u.aln",
      prefxCStr,
      alnST->refStartUI,
      alnST->refEndUI,
      alnST->queryStartUI,
      alnST->queryEndUI
    ); // Make the new file name
  
    outFILE = fopen(fileNameCStr, "w");
  
    if(outFILE == 0)
    { // If I could not open the file
      free(refAlnCStr);
      free(queryAlnCStr);
      freeAlnST(alnST, 1); // No longer need
  
      return 1;
    } // If I could not open the file
  
    // Print the alignment
    printAln(
      outFILE,
      queryST->idCStr,
      queryAlnCStr,
      refAlnCStr,
      refST->idCStr,
      alnMtxST->bestScoreST.scoreL,
      settings->lineWrapUS,
      alnST
    );
  
    /******************************************************\
    * Fun-03 Sec-05 Sub-05:
    *  - Free uneeded variables (clean up)
    \******************************************************/
  
    free(refAlnCStr);
    refAlnCStr = 0;
  
    free(queryAlnCStr);
    queryAlnCStr = 0;
  
    freeAlnST(alnST, 1); // No longer need
    alnST = 0;
  
    fclose(outFILE);
    ++scoreST;  // Move to the next entry
  } // For all reference bases in the alignment

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
   ' Fun-04 TOC: updateDirScoreWaterSingle
   '  - Picks the best score and direction for the current
   '    base pairs being compared in a Waterman Smith
   '    alignment
   '  o fun-04 sec-1: Matches->insertions->deletions
   '  o fun-04 sec-2: Matches->deletions->insertions
   '  o fun-04 sec-3: Insertions->matches->deletions
   '  o fun-04 sec-4: Deletions->matches->insertions
   '  o fun-04 sec-5: Insertions->deletions->matches
   '  o fun-04 sec-6: Deletions->insertions->matches
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-1: Matches and then insertions
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
               ^ Fun-04 Sec-2: Matches->dels->ins
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
       ^ Fun-04 Sec-3: Insertions->matches->deletions
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
               ^ Fun-04 Sec-4: Dels->matches->ins
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
       ^ Fun-04 Sec-5: Insertions->deletions->matches
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
               ^ Fun-04 Sec-6: Dels->insertions->matches
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
