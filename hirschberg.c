/*#########################################################
# Name: hirschberg
# Use:
#  - Holds functions for doing a hirschberg global
#    alignment
# Libraries:
#  - neelde.h
#  - water.h
#  o "generalAlnFun.h"
#  o "alnStruct.h"
#  o "alnMatrixStruct.h"
#  o "twoBitArrays.h"
#  o "scoresST.h"
#  o "seqStruct.h"
#  o "alnSetStruct.h"
#  o "alnSeqDefaults.h"
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>  // by alnSetStructure.h
#  o <string.h>
#########################################################*/

#include "hirschberg.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 Hirschberg:
'    - Sets up for and calls the recursvie function to
'      run a Hirschberg alignment
'  o fun-02 HirschbergFun:
'    - Does the recursive part of a Hirschberg alignment
'  o fun-03 scoreForwardHirsch:
'    - Does a single round of scoring for a hirschberg
'      alignment (forward direction)
'  o fun-04 scoreReverseHirsch:
'    - Does a single round of scoring for a hirschberg
'      alignment (reverse direction)
'  o fun-05 positionSingleRefBase:
'    - Align a single reference base to a query sequence
'  o fun-06 twoBitAlnAryToAlnST:
'    - Converts a two bit array with an alignment to an
'      alnStruct structure for printing
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o A alignment structure with the alignment.
|    o 0 For memory errors
\--------------------------------------------------------*/
struct alnStruct * Hirschberg(
  struct seqStruct *refST, // Reference sequence to align
  struct seqStruct *qryST, // Qeury sequence to align
    // For refST and qryST, use seqStruct->offsetUL to set
    // the starting point for the alignmnet and
    // seqStruct->endAlnUL to set the ending point
  struct alnSet *settings  // Settings to use for alignment
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Hirschberg
   '  - Sets up for and calls the recursvie function to
   '    run a Hirschberg alignment
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Memory allocation (set up for Hirschberg)
   '  o fun-01 sec-03:
   '    - Run the hirschberg alignment
   '  o fun-01 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   unsigned long lenRefUL =
     refST->endAlnUL - refST->offsetUL + 1;
     // +1 to convert to index 1 (values are index 0)
   unsigned long lenQryUL =
     qryST->endAlnUL - qryST->offsetUL + 1;
     // + 1 to convert to index 1 (values are index 0)

   long *forwardScoreRowL = 0;
   long *reverseScoreRowL = 0;

   struct twoBitAry *refTwoBitST = 0;// Holds ref alignment
   struct twoBitAry *qryTwoBitST = 0;//Holds queryAlignment
   //struct alnStruct *alignmentST = 0; // returned alignment

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-02:
   ^   - Memory allocation (set up for Hirschberg)
   ^   o fun-01 sec-02 sub-01:
   ^     - Initalize the ouput alignment structure 
   ^   o fun-01 sec-02 sub-02:
   ^     - Initalize the scoring rows
   ^   o fun-01 sec-02 sub-03:
   ^     - Initalize the direction rows
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Initalize the ouput alignment structure 
   \******************************************************/

   refTwoBitST = makeTwoBitArry(lenRefUL, 0);
    // + 2 for index 0

   if(refTwoBitST == 0) return 0;

   // Makign sure that the query alignment array is at
   // least as long as tle reference. This is so I can
   // use in the reverse scoring and keepl alnSeq thread
   // safe
   if(lenQryUL >= lenRefUL)
     qryTwoBitST = makeTwoBitArry(lenQryUL, 0);
   else
     qryTwoBitST = makeTwoBitArry(lenRefUL, 0);

   if(qryTwoBitST == 0)
   { // If had a memroy allocation error
     freeTwoBitAry(refTwoBitST, 0, 0);
     return 0;
   } // If had a memroy allocation error

   /******************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Initalize the scoring rows
   \******************************************************/

   // I am using full length arrays to make the later
   // steps eaiser. This takes more memory, but makes life
   // nicer

   forwardScoreRowL = malloc(sizeof(long) * lenRefUL);

   if(forwardScoreRowL == 0)
   { // If had a memory allocatoin error
     freeTwoBitAry(refTwoBitST, 0, 0);
     freeTwoBitAry(qryTwoBitST, 0, 0);
     return 0;
   } // If had a memory allocatoin error

   reverseScoreRowL = malloc(sizeof(long) * lenRefUL);

   if(reverseScoreRowL == 0)
   { // If had a memory allocatoin error
     freeTwoBitAry(refTwoBitST, 0, 0);
     freeTwoBitAry(qryTwoBitST, 0, 0);
     free(forwardScoreRowL);
     return 0;
   } // If had a memory allocatoin error

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-03:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   HirschbergFun(
     refST->seqCStr + refST->offsetUL,
     0,                // Sendingi in the starting base
     lenRefUL,         // Length of ref region to align
     qryST->seqCStr + qryST->offsetUL,
     0,                // Sending 1st query base to align
     lenQryUL,         // length of query target region
     forwardScoreRowL, // For scoring
     reverseScoreRowL, // For scoring
     refTwoBitST,      // HOlds the reference alignment
     qryTwoBitST,      // Holds the query alignment
     settings         // Settings for the alignment
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-04:
   ^    - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(forwardScoreRowL);
   free(reverseScoreRowL);

   // Free the variables that are no longer used
   printTwoBitAln(
     stdout,
     refST->seqCStr,
     qryST->seqCStr,
     refTwoBitST,
     qryTwoBitST,
     59
   );

   freeTwoBitAry(refTwoBitST, 0, 0);
   freeTwoBitAry(qryTwoBitST, 0, 0);

   return 0;
} // Hirschberg

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the output alignment
\--------------------------------------------------------*/
void HirschbergFun(
  char *refSeqCStr,          // Reference sequence
  unsigned long refStartUL,// index 0 Starting point on ref
  unsigned long refLenUL,//index 1 length of region toAlign

  char *qrySeqCStr,          // Query sequence
  unsigned long qryStartUL,//index 0 Starting query base
  unsigned long qryLenUL,//index 1 Length of target

  long *forwardScoreRowL,   // Holds final forward row
  long *reverseScoreRowL,   // For finding reverse scores
    // both the forward and reverse scoring rows must be
    // the size of the full length reference.
  struct twoBitAry *refAlnST,   // For gap extension
  struct twoBitAry *qryAlnST,  // Holds alignment codes
  struct alnSet *settings  // Settings to use for alignment
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: HirschbergFun
   '  - Does the recursive part of a Hirschberg alignment
   '  o fun-02 sec-01:
   '    - Variable declerations
   '  o fun-02 sec-02:
   '    - Check if on a leaf (final part of alignment
   '  o fun-02 sec-03:
   '    - Get scores
   '  o fun-02 sec-04:
   '    - Find the midpoint
   '  o fun-02 sec-05:
   '    - Run the next hirschberg alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint8_t bitC = 0;
   long forwardIndelColL = 0;
   long reverseIndelColL = 0;
   unsigned long midPointUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Check if on a leaf (final part of alignment
   ^  o fun-02 sec-02 sub-01:
   ^    - Handle cases were I have just insertions
   ^  o fun-02 sec-02 sub-02:
   ^    - Handle cases were I have just deletions
   ^  o fun-02 sec-02 sub-03:
   ^    - Handle cases were I have to align last ref base
   ^  o fun-02 sec-02 sub-04:
   ^  - Handle cases were I have to align last query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-01:
   *  - Handle cases were I have just insertions
   \******************************************************/

   if(refLenUL == 0)
   { // If all remaing bases are query insertions
     moveXElmFromStart(qryAlnST, qryStartUL);

     for(
       unsigned long numInsUI = 0;
       numInsUI < qryLenUL;
       ++numInsUI
     ){ // Loop: fill in the insertions
       changeTwoBitElm(qryAlnST, defGapFlag);
       twoBitAryMoveToNextElm(qryAlnST);
     } // Loop: fill in the insertions

     return; // Nothing else to do
   } // If all remaing bases are query insertions

   /******************************************************\
   * Fun-02 Sec-02 Sub-02:
   *  - Handle cases were I have just deletions
   \******************************************************/

   if(qryLenUL == 0)
   { // If all remaing bases are query deletions
     moveXElmFromStart(refAlnST, refStartUL);

     for(
       unsigned long numInsUI = 0;
       numInsUI < refLenUL;
       ++numInsUI
     ){ // Loop: fill in the insertions
       changeTwoBitElm(refAlnST, defGapFlag);
       twoBitAryMoveToNextElm(refAlnST);
     } // Loop: fill in the insertions

     return; // Nothing else to do
   } // If all remaing bases are query deletions

   /******************************************************\
   * Fun-02 Sec-02 Sub-03:
   *  - Handle cases were I have to align last ref base
   \******************************************************/

   if(refLenUL == 1)
   { // If I have to align the last reference base
     moveXElmFromStart(qryAlnST, qryStartUL);
     moveXElmFromStart(refAlnST, refStartUL);

     if(qryLenUL == 0)
     { // If bases are aligned (one reference & one query)
        changeTwoBitElm(refAlnST, defGapFlag);
        return; // Finished
     } // If bases are aligned (one reference & one query)

     if(qryLenUL == 1)
     { // If bases are aligned (one reference & one query)
       qrySeqCStr += qryStartUL;
       refSeqCStr += refStartUL;

       if(checkIfBasesMatch(qrySeqCStr, refSeqCStr) == 1)
         bitC = defMatchFlag;

        else bitC = defSnpFlag;

        changeTwoBitElm(qryAlnST, bitC);
        changeTwoBitElm(refAlnST, bitC);
        return; // Finished
     } // If bases are aligned (one reference & one query)

     positionSingleBase(
       *(refSeqCStr + refStartUL),// ref base
       refStartUL,                // Position of ref base
       qrySeqCStr,              // first base of query
       qryStartUL,              // positoin of query
       qryLenUL,                // Length of the query
       refAlnST,                // Array to hold alignment
       qryAlnST,                // Array to hold alignment
       settings                 // Has Scoring variables
     );

     return; // This base is now aligned
   } // If I have to align the last reference base

   /******************************************************\
   * Fun-02 Sec-02 Sub-04:
   *  - Handle cases were I have to align last query base
   \******************************************************/

   if(qryLenUL == 1)
   { // If I have to align the last query base
     moveXElmFromStart(qryAlnST, qryStartUL);
     moveXElmFromStart(refAlnST, refStartUL);

     if(refLenUL == 0)
     { // If bases are aligned (one reference & one query)
        changeTwoBitElm(qryAlnST, defGapFlag);
        return; // Finished
     } // If bases are aligned (one reference & one query)

     if(refLenUL == 1)
     { // If bases are aligned (one reference & one query)
       qrySeqCStr += qryStartUL;
       refSeqCStr += refStartUL;

       if(checkIfBasesMatch(qrySeqCStr, refSeqCStr) == 1)
         bitC = defMatchFlag;

        else bitC = defSnpFlag;

        changeTwoBitElm(qryAlnST, bitC);
        changeTwoBitElm(refAlnST, bitC);
        return; // Finished
     } // If bases are aligned (one reference & one query)

     positionSingleBase(
       *(qrySeqCStr + qryStartUL),// ref base
       qryStartUL,                // Position of ref base
       refSeqCStr,              // first base of reference
       refStartUL,              // positoin of query
       refLenUL,                // Length of the query
       qryAlnST,                // Array to hold alignment
       refAlnST,                // Array to hold alignment
       settings                 // Has Scoring variables
     );

     return; // Finshed aligning this query base
   } // If I have to align the last query base

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-03:
   ^  - Get scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // For multithreading I would need need a third
    // direction row. Otherwise I will have thread racing
    // for the directional matrix
    forwardIndelColL = 
      scoreForwardHirsch(
        refSeqCStr,       // Entire reference sequence
        refStartUL,       // Starting base of ref target
        refLenUL,         // length of ref target region
        qrySeqCStr,       // Query seq with coordinates
        qryStartUL,     // Starting base of query target
        qryLenUL / 2,   // Length of query target region
        forwardScoreRowL, // Array of scores to fill
        refAlnST,         // direction row for gap extend
        settings         // setttings to use
    ); // Get the scores for the forward direction

    reverseIndelColL = 
      scoreReverseHirsch(
        refSeqCStr,         // Entire reference sequence
        refStartUL,         // Starting base of ref target
        refLenUL,           // length of ref target region
        qrySeqCStr,         // Query seq with coordinates
        qryStartUL + (qryLenUL / 2),//Start of query target
        qryLenUL - (qryLenUL / 2), // New query length
        reverseScoreRowL,   // Array of scores to fill
        refAlnST,           // direction row for gap extend
        settings            // setttings to use
      ); // Get the scores for the reverse direction
      // I can get away with queryLen/2 here, because 
      // queryLen is index 1 and the function takes in
      // an lenth 1 argument

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-04:
   ^   - Find the midpoint
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *(forwardScoreRowL + refStartUL + refLenUL - 1) +=
     reverseIndelColL;

   midPointUL = refLenUL;

   for(
     unsigned long baseUL = 0;
     baseUL < refLenUL - 1;
     ++baseUL
   ) { // Loop; add up all scores
     *(forwardScoreRowL + refStartUL + baseUL) += 
       *(reverseScoreRowL + refStartUL + baseUL + 1);
       // The reverse row is already reversed

     if(
       *(forwardScoreRowL + refStartUL + baseUL) >
       *(forwardScoreRowL + refStartUL + midPointUL -1)
     ) midPointUL = baseUL + 1;
   }// Loop; add up all scores

   *(reverseScoreRowL + refStartUL) += forwardIndelColL;

   if(
     *(reverseScoreRowL + refStartUL) >
     *(forwardScoreRowL + refStartUL + midPointUL - 1)
   ) midPointUL = 0;


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-05:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   HirschbergFun(
     refSeqCStr,      // Full reference sequence
     refStartUL,      // Full reference sequence
     midPointUL,
     qrySeqCStr,      // Full queyr sequence
     qryStartUL,    // Start of queyr target region
     qryLenUL / 2,    // Length of query target region
     forwardScoreRowL,  // For scoring
     reverseScoreRowL, // Has last line of scores
     refAlnST,        // direction row for gap extend
     qryAlnST,        // Holds the alignment codes
     settings         // Settings for the alignment
   );

   HirschbergFun(
     refSeqCStr,                 // Full reference sequence
     refStartUL + midPointUL,    // New reference start
     refLenUL - midPointUL,        // New reference end
     qrySeqCStr,                  // Full query sequence
     qryStartUL + (qryLenUL / 2),
     qryLenUL - (qryLenUL / 2),     // New query length
     forwardScoreRowL,  // For scoring
     reverseScoreRowL, // Has last line of scores
     refAlnST,        // direction row for gap extend
     qryAlnST,        // Holds the alignment codes
     settings         // Settings for the alignment
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-0?:
   ^    - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return;
} // HirschbergFun

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      Needleman Wunsch / Smith Waterman alignment
|    o dirRowSt to hold the last row of directions in a
|      Needleman Wunsch / Smith Waterman alignment
\--------------------------------------------------------*/
long scoreForwardHirsch(
  char *refSeqCStr,          // Reference sequence
  unsigned long refStartUL,  // index 0 starting ref base
  unsigned long refLenUL,    // index 1 Length of target

  char *qrySeqCStr,          // Query sequence
  unsigned long qryStartUL,// Index 0 Starting query base
  unsigned long qryLenUL,  // index 1 length of target

  long *scoreRowPtrL,        // Array of scores to fill
  struct twoBitAry *dirRowST,//direction row for gap extend
  struct alnSet *settings    // setttings to use
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: scoreForwardHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (forward direction)
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-03 sec-03:
   '    - Score till on the last row
   '  o fun-03 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long *scoreOnPtrL = scoreRowPtrL + refStartUL;
   long indelColL = 0;
   long insScoreL = 0;
   long delScoreL = 0;
   long matchScoreL = 0;
   long nextMatchScoreL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Set up the first row (indel row) of scores
   ^  o fun-03 sec-02 sub-01:
   ^    - Set up the first two elements (no gat extend)
   ^  o fun-03 sec-02 sub-02:
   ^    - Set up remaing elements (indel) (uses gap extend)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - Set up the first two elements (no gat extend)
   \******************************************************/

   moveXElmFromStart(dirRowST, refStartUL);

   indelColL = 0;

   *scoreOnPtrL = settings->gapStartPenaltyI;
   changeTwoBitElm(dirRowST, defMoveLeft);

   ++scoreOnPtrL;
   twoBitAryMoveToNextElm(dirRowST);

   /******************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Set up remaing elements (indels) (uses gap extend)
   \******************************************************/

   for(
     char *refBaseCStr = refSeqCStr + refStartUL + 1;
     refBaseCStr < refSeqCStr + refStartUL + refLenUL;
     ++refBaseCStr
   ){ // Loop:Set the initial blank scores
     *scoreOnPtrL =
       *(scoreOnPtrL-1) +settings->gapExtendPenaltyI;

     changeTwoBitElm(dirRowST, defMoveLeft);

     ++scoreOnPtrL;
     twoBitAryMoveToNextElm(dirRowST);
   } // Loop:Set the initial blank scores

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Score till on the last row
   ^  o fun-03 sec-03 sub-01:
   ^    - Set up the first element (indel, has gap open)
   ^  o fun-03 sec-03 sub-02:
   ^    - Find the scores for the next row
   ^  o fun-03 sec-03 sub-03:
   ^    - Select best score for direction
   ^  o fun-03 sec-03 sub-04:
   ^    - Set up for scoring next row (handle indel col)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-01:
   *  - Set up the first element (indel, has gap open)
   \******************************************************/

   // Move back to the start of row
   moveXElmFromStart(dirRowST, refStartUL);
   scoreOnPtrL = scoreRowPtrL + refStartUL;

   // Fine the first match/snp score (first ref base)
   nextMatchScoreL =
     getBasePairScore(
       qrySeqCStr + qryStartUL,   // first query base
       refSeqCStr + refStartUL,     // first ref base
       settings                     // Has score matrix
   ); // Ge the score for the enxt base

   nextMatchScoreL += indelColL;
     // This is here so I can overwrite the array with the
     // new scores

   indelColL += settings->gapStartPenaltyI;

   // Find the first insertion and deletion scores
   insScoreL =  *scoreOnPtrL + settings->gapExtendPenaltyI;
   delScoreL = indelColL + settings->gapExtendPenaltyI;

   /******************************************************\
   * Fun-03 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   for(
     char *qryBaseCStr = qrySeqCStr + qryStartUL;
     qryBaseCStr < qrySeqCStr + qryStartUL + qryLenUL;
     ++qryBaseCStr
   ){ // Loop: score all query bases (rows)

     for(
       char *refBaseCStr = refSeqCStr + refStartUL;
       refBaseCStr < refSeqCStr + refStartUL + refLenUL -1;
       ++refBaseCStr
     ){ // Loop:score all query bases (columns)

       // Get the score for the next match. This allows me
       // to overwite the last diagnol and thus, use only
       // one row of scoring

       matchScoreL = nextMatchScoreL;

       nextMatchScoreL =
         getBasePairScore(
           qryBaseCStr,     // Current query base
           refBaseCStr + 1, // next reference base
           settings         // Has score matrix
       ); // Ge the score for the enxt base
         // At worst case refBaseCStr will be '\0'

       nextMatchScoreL += *scoreOnPtrL;

       /**************************************************\
       * Fun-03 Sec-03 Sub-03:
       *  - Select best score for direction
       \**************************************************/

       updateDirAndScoreNeedle(
         dirRowST,        // Direction matrix
         settings,        // has direction preference
         &insScoreL,      // Score for an insertion
         &matchScoreL,    // Score for an deletion
         &delScoreL,      // Score for an match/snp
         scoreOnPtrL        // Score position to update
       ); // Update the score using a waterman

       // The deletion scores are based on the found base,
       // So I can find the next score before moving
       // Get the deletion score
       delScoreL =
         getIndelScore(dirRowST, settings, scoreOnPtrL);
         // This is fun-03 in generalAlnFun.h

       // Move to the next element
       twoBitAryMoveToNextElm(dirRowST);
       ++scoreOnPtrL;

       // Finding indel scores at end, so that I can keep
       // the indel column in a separate variable
       // Get the insertion score
       insScoreL =
         getIndelScore(dirRowST, settings, scoreOnPtrL);
     } // Loop:score all query bases (columns)

     /****************************************************\
     * Fun-03 Sec-03 Sub-04:
     *  - Set up for scoring next row (handle indel col)
     \****************************************************/

     // Update the final score in the row
     updateDirAndScoreNeedle(
       dirRowST,        // Direction matrix
       settings,        // has direction preference
       &insScoreL,      // Score for an insertion
       &nextMatchScoreL,    // Score for an deletion
       &delScoreL,      // Score for an match/snp
       scoreOnPtrL        // Score position to update
     ); // Update the score using a waterman

     // Move back to the start of row
     moveXElmFromStart(dirRowST, refStartUL);
     scoreOnPtrL = scoreRowPtrL + refStartUL;

     /// I need to refind the insertion and deletion scores
     insScoreL =
       getIndelScore(dirRowST, settings, scoreOnPtrL);

     // Fine the first match/snp score (first ref base)
     nextMatchScoreL =
       getBasePairScore(
         qryBaseCStr + 1,           // Next query base
         refSeqCStr + refStartUL,   // first base
         settings                   // Has score matrix
     ); // Ge the score for the enxt base

     nextMatchScoreL += indelColL;
       // This is here so I can overwrite the array with the
       // new scores

     indelColL += settings->gapExtendPenaltyI;
     delScoreL = indelColL+settings->gapExtendPenaltyI;
   } // Loop: score all query bases (rows)
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Clean up
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return indelColL;
} // scoreForwardHirsch

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      backwards Needleman Wunsch /Smith Waterman alignment
|    o dirRowST to hold the last row of directions in a
|      backwards Needleman Wunsch /Smith Waterman alignment
\--------------------------------------------------------*/
long scoreReverseHirsch(
  char *refSeqCStr,          // Reference sequence
  unsigned long refStartUL,  // index 0 starting ref base
  unsigned long refLenUL,    // index 1 Length of target

  char *qrySeqCStr,          // Query sequence
  unsigned long qryStartUL,// Index 0 Starting query base
  unsigned long qryLenUL,  // index 1 length of target

  long *scoreRowPtrL,       // Array of scores to fill
     // This needs to be as long as the full length
     // reference sequence
  struct twoBitAry *dirRowST,//direction row for gap extend
  struct alnSet *settings   // setttings to use
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: scoreReverseHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (reverse direction)
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-04 sec-03:
   '    - Score till on the last row
   '  o fun-04 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long *scoreOnPtrL = 0;
     // refLenUL is index 1, so  need -1 to get index 0

   long insScoreL = 0;
   long delScoreL = 0;
   long matchScoreL = 0;
   long nextMatchScoreL = 0; // Score for the next base
   long indelColL = 0;       // Holds indel column values

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Set up the first row (indel row) of scores
   ^  o fun-04 sec-02 sub-01:
   ^    - Set up the first two elements (no gat extend)
   ^  o fun-04 sec-02 sub-02:
   ^    - Set up remaing elements (indel) (uses gap extend)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-02 Sub-01:
   *  - Set up the first two elements (no gat extend)
   \******************************************************/

   moveXElmFromStart(dirRowST, refStartUL + refLenUL - 1);
   scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL - 1;
     // - 1 to account for refLenUL being index 1


   indelColL = 0;

   *scoreOnPtrL = settings->gapStartPenaltyI;
   changeTwoBitElm(dirRowST, defMoveLeft);

   --scoreOnPtrL;
   twoBitAryMoveBackOneElm(dirRowST);

   /******************************************************\
   * Fun-04 Sec-02 Sub-02:
   *  - Set up remaing elements (indels) (uses gap extend)
   \******************************************************/

   // Loop from the second to last base till the start,
   // Starting at second to last, because already did the
   // first base
   for(
     char *refBaseCStr =
       refSeqCStr + refStartUL + refLenUL - 2;
     refBaseCStr > refSeqCStr + refStartUL - 1; // index 0
     --refBaseCStr
   ){ // Loop:Set the initial blank scores
     *scoreOnPtrL =
       *(scoreOnPtrL+1) +settings->gapExtendPenaltyI;

     changeTwoBitElm(dirRowST, defMoveLeft);

     --scoreOnPtrL;
     twoBitAryMoveBackOneElm(dirRowST);
   } // Loop:Set the initial blank scores

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Score till on the last row
   ^  o fun-04 sec-03 sub-01:
   ^    - Set up the first element (indel, has gap open)
   ^  o fun-04 sec-03 sub-02:
   ^    - Find the scores for the next row
   ^  o fun-04 sec-03 sub-03:
   ^    - Select best score for direction
   ^  o fun-04 sec-03 sub-04:
   ^    - Set up for scoring next row (handle indel col)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-03 Sub-01:
   *  - Set up the first element (indel, has gap open)
   \******************************************************/

   // Move back to the start of row (refLenUL is index 1)
   moveXElmFromStart(dirRowST, refStartUL + refLenUL - 1);
   scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL - 1;

   // Fine the first match/snp score (first ref base)
   nextMatchScoreL =
     getBasePairScore(
       qrySeqCStr + qryStartUL + qryLenUL - 1,
       refSeqCStr + refStartUL + refLenUL - 1,
       settings                         // Has score matrix
   ); // Ge the score for the enxt base

   nextMatchScoreL += indelColL;
     // This is here so I can overwrite the array with the
     // new scores

   // Find the first insertion and deletion scores
   indelColL += settings->gapStartPenaltyI;
   insScoreL = indelColL + settings->gapExtendPenaltyI;
   delScoreL = indelColL + settings->gapExtendPenaltyI;

   /******************************************************\
   * Fun-04 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   for(
     char *qryBaseCStr =
       qrySeqCStr + qryStartUL + qryLenUL - 1;
     qryBaseCStr > qrySeqCStr + qryStartUL - 1;
     --qryBaseCStr
   ){ // Loop: score all query bases (rows)

     for(
       char *refBaseCStr =
         refSeqCStr + refStartUL + refLenUL - 1;
       refBaseCStr > refSeqCStr + refStartUL;
       --refBaseCStr
     ){ // Loop:score all query bases (columns)

       // Get the score for the next match. This allows me
       // to overwite the last diagnol and thus, use only
       // one row of scoring

       matchScoreL = nextMatchScoreL;

       nextMatchScoreL =
         getBasePairScore(
           qryBaseCStr,     // Current query base
           refBaseCStr - 1, // next reference base
           settings         // Has score matrix
       ); // Ge the score for the enxt base
         // At worst case refBaseCStr will be '\0'

       nextMatchScoreL += *scoreOnPtrL;

       /**************************************************\
       * Fun-04 Sec-03 Sub-03:
       *  - Select best score for direction
       \**************************************************/

       updateDirAndScoreNeedle(
          dirRowST,        // Direction matrix
          settings,        // has direction preference
          &insScoreL,      // Score for an insertion
          &matchScoreL,    // Score for an deletion
          &delScoreL,      // Score for an match/snp
          scoreOnPtrL        // Score position to update
        ); // Update the score using a waterman

       // The deletion scores are based on the found base,
       // So I can find the next score before moving
       // Get the deletion score
       delScoreL =
         getIndelScore(dirRowST, settings, scoreOnPtrL);
         // This is fun-03 in generalAlnFun.h

       twoBitAryMoveBackOneElm(dirRowST);
       --scoreOnPtrL;

       // Finding indel scores at end, so that I can keep
       // the indel column in a separate variable
       // Get the insertion score
       insScoreL =
         getIndelScore(dirRowST, settings, scoreOnPtrL);
     } // Loop:score all query bases (columns)

     /****************************************************\
     * Fun-04 Sec-03 Sub-04:
     *  - Set up for scoring next row (handle indel col)
     \****************************************************/

     // Update the final cell in the row
     updateDirAndScoreNeedle(
        dirRowST,        // Direction matrix
        settings,        // has direction preference
        &insScoreL,      // Score for an insertion
        &nextMatchScoreL,    // Score for an deletion
        &delScoreL,      // Score for an match/snp
        scoreOnPtrL        // Score position to update
      ); // Update the score using a waterman

     // Move back to the start of row
     moveXElmFromStart(dirRowST, refStartUL + refLenUL -1);
     scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL -1;

     /// I need to refind the insertion and deletion scores
     insScoreL =
       getIndelScore(dirRowST, settings, scoreOnPtrL);

     // Fine the first match/snp score (first ref base)
     nextMatchScoreL =
       getBasePairScore(
         qryBaseCStr - 1,           // Next query base
         refSeqCStr + refStartUL + refLenUL - 1,
         settings         // Has score matrix
     ); // Ge the score for the enxt base

     nextMatchScoreL += indelColL;
       // This is here so I can overwrite the array with the
       // new scores

     indelColL += settings->gapExtendPenaltyI;
     delScoreL = indelColL+settings->gapExtendPenaltyI;
   } // Loop: score all query bases (rows)
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-04:
  ^  - Clean up
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return indelColL;
} // scoreReverseHirsch

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the alignment for the single
|      base aligned to the sequence
\--------------------------------------------------------*/
void positionSingleBase(
  char baseC,             // Single base to align to a seq
  unsigned long baseIndexUL, // Index base is at
  char *seqCStr,            // Sequence to position base on
  unsigned long startOfSeqUL,
    // Index 0 of first base to align bascC to in seqCStr
  unsigned long lenSeqUL,
    //index 1; Length of the aligned region in seqCStr
  struct twoBitAry *baseCAlnST,
    // Two bit alingment array for the sequence having
    // baseC
  struct twoBitAry *seqAlnST,    // Array to hold alignment
    // Two bit alignment array for the sequence aliging
    // baseC to
  struct alnSet *settings        // setttings to use
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBase
   '  - Align a single base to a sequence
   '  o fun-05 sec-01:
   '    - Variable declerations
   '  o fun-05 sec-02:
   '    - Find the reference bases position on the query
   '  o fun-05 sec-03:
   '    - Fill in insertions and reference base position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long insScoreL = 0;      // The score of an insertion
   long delScoreL = 0;      // the score of a deleton
   long matchScoreL = 0;    // The score of a match
   long curScoreL = 0;      // The score of the last row

   char *seqBaseCStr = seqCStr + startOfSeqUL;

   char *endSeqCStr = seqCStr +startOfSeqUL +lenSeqUL-1;

   char *lastMatchCStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-02:
   ^  - Find the reference bases position on the query
   ^  o fun-05 sec-02 sub-01:
   ^    - Find the first scores for the loop
   ^  o fun-05 sec-02 sub-02:
   ^    - Find the remaing scores
   ^  o fun-05 sec-02 sub-03:
   ^     - Figure out which of hte final scores to keep
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-05 Sec-02 Sub-01:
   *  - Find the first scores for the loop
   \******************************************************/

   moveXElmFromStart(baseCAlnST, baseIndexUL);
   moveXElmFromStart(seqAlnST, startOfSeqUL);

   matchScoreL =
       getBasePairScore(seqBaseCStr, &baseC, settings);

   insScoreL = settings->gapStartPenaltyI;
   delScoreL = insScoreL;

   ++seqBaseCStr;

   /******************************************************\
   * Fun-05 Sec-02 Sub-02:
   *  - Find the remaing scores
   \******************************************************/

 
   do { // While I have bases to compare
     // Set your score and direction
     updateDirAndScoreNeedle(
       seqAlnST,        // Direction matrix
       settings,        // has direction preference
       &insScoreL,      // Score for an insertion
       &matchScoreL,    // Score for an deletion
       &delScoreL,      // Score for an match/snp
       &curScoreL        // Score position to update
     ); // Update the score using a waterman

     matchScoreL =
         insScoreL
       + getBasePairScore(seqBaseCStr,&baseC,settings);

     insScoreL += settings->gapExtendPenaltyI;

     switch(getTwoBitAryElm(seqAlnST))
     { // Switch: Check if keeping the score
       case defMoveDiagnol:
         lastMatchCStr = seqBaseCStr - 1;
         delScoreL = curScoreL+settings->gapStartPenaltyI;
         break;

       case defMoveLeft:
       case defMoveUp:
         delScoreL = curScoreL+settings->gapExtendPenaltyI;
         break;
       case defMoveStop: break; // Never fires
     } // Switch: Check if keeping the score

     ++seqBaseCStr;
   } while(seqBaseCStr <= endSeqCStr);
   // While I have bases to compare

   /******************************************************\
   * Fun-05 Sec-02 Sub-03:
   *  - Figure out which of hte final scores to keep
   \******************************************************/

   // Set your score and direction
   updateDirAndScoreNeedle(
     seqAlnST,        // Direction matrix
     settings,        // has direction preference
     &insScoreL,      // Score for an insertion
     &matchScoreL,    // Score for an deletion
     &delScoreL,      // Score for an match/snp
     &curScoreL        // Score position to update
   ); // Update the score using a waterman

   switch(getTwoBitAryElm(seqAlnST))
   { // Switch: Check if keeping the score
       case defMoveDiagnol:
         lastMatchCStr = seqBaseCStr - 1;
         break;

       case defMoveLeft: break;
       case defMoveUp: break;
       case defMoveStop: break; // Never fires
   } // Switch: Check if keeping the score

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-03:
   ^  - Fill in the insertions and reference base position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   // A series of deletions and insertions are prefered
   // over matches and smps. In this case put the base
   // at the first base. There is no good position
   if(lastMatchCStr == 0)
     lastMatchCStr = seqCStr + startOfSeqUL;

   seqBaseCStr = seqCStr + startOfSeqUL;
   // No need to change query position since previouis
   // loop only usedo one direction position

   // Add in the insertions at the start
   while(seqBaseCStr < lastMatchCStr)
   { // While I have insertions to fill
     changeTwoBitElm(seqAlnST, defGapFlag);
     twoBitAryMoveToNextElm(seqAlnST);
     ++seqBaseCStr;
   } // While I have insertions to fill
   
   // Add in the positoin of the base
   if(checkIfBasesMatch(seqBaseCStr, &baseC))
   { // IF the bases matched
      changeTwoBitElm(baseCAlnST, defMatchFlag);
      changeTwoBitElm(seqAlnST, defMatchFlag);
   } // IF the bases matched

    else
    { // Else this was a SNP
      changeTwoBitElm(baseCAlnST, defSnpFlag);
      changeTwoBitElm(seqAlnST, defSnpFlag);
    } // Else this was a SNP

   twoBitAryMoveToNextElm(seqAlnST);
   ++seqBaseCStr;

   // Finish adding in the insertions at the end
   while(seqBaseCStr <= endSeqCStr)
   { // While I have insertions to fill
     changeTwoBitElm(seqAlnST, defGapFlag);
     twoBitAryMoveToNextElm(seqAlnST);
     ++seqBaseCStr;
   } // While I have insertions to fill

   return;
} // positionSingleRefBase

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o 0: for success
|    o 1: No valid output file (alnFILE)
|    o 64: for memroy error
|  - Prints:
|    o Alginmetns to alnFILE
\--------------------------------------------------------*/
char printTwoBitAln(
  FILE *alnFILE,               // File to save alignment to
  char *refSeqCStr,
  char *qrySeqCStr,
  struct twoBitAry *refTwoBitST,
    // Two bir array with the referenc alignmetn
  struct twoBitAry *qryTwoBitST,
    // Two bir array with the query alignmetn
  unsigned long lnWrapUL
    // How many characters to pritn per line. Do
    // length(query) + length(ref) to have no line wrap
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: printTwoBitAln
   '  - Converts a two bit array with an alignment to an
   '    alnStruct structure for printing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *refAlnCStr = 0;
   char *qryAlnCStr = 0;
   char *eqxAlnCStr = 0;

   char *tmpRefCStr = 0;
   char *tmpQryCStr = 0;
   char *tmpEqxCStr = 0;

   unsigned long lenBuffUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-02:
   ^  - Allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   moveXElmFromStart(refTwoBitST, 0);
   moveXElmFromStart(qryTwoBitST, 0);

   if(alnFILE == 0) return 1; // Invalid file

   refAlnCStr = malloc(sizeof(char) * (lnWrapUL + 2));

   if(refAlnCStr == 0) return 64;

   qryAlnCStr = malloc(sizeof(char) * (lnWrapUL + 2));

   if(qryAlnCStr == 0)
   { // IF had a memory allocation error
     free(refAlnCStr);
     refAlnCStr = 0;
     return 64;
   } // IF had a memory allocation error

   eqxAlnCStr = malloc(sizeof(char) * (lnWrapUL + 2));

   if(eqxAlnCStr == 0)
   { // IF had a memory allocation error
     free(refAlnCStr);
     refAlnCStr = 0;

     free(eqxAlnCStr);
     eqxAlnCStr = 0;

     return 64;
   } // IF had a memory allocation error

   *(refAlnCStr + lnWrapUL) = '\n';
   *(qryAlnCStr + lnWrapUL) = '\n';
   *(eqxAlnCStr + lnWrapUL) = '\n';

   *(refAlnCStr + lnWrapUL + 1) = '\0';
   *(qryAlnCStr + lnWrapUL + 1) = '\0';
   *(eqxAlnCStr + lnWrapUL + 1) = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-03:
   ^  - Finish printing the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tmpRefCStr = refAlnCStr;
   tmpQryCStr = qryAlnCStr;
   tmpEqxCStr = eqxAlnCStr;

    while(*refSeqCStr != '\0' && *qrySeqCStr != '\0')
    { // While I have gaps at the start

      // Check if I have any insertions to handle
      if(getTwoBitAryElm(qryTwoBitST) == defGapFlag)
      { // If have an insertion
        *tmpRefCStr = '-';
        *tmpQryCStr = *qrySeqCStr;
        *tmpEqxCStr = 'I';

        ++qrySeqCStr;
        twoBitAryMoveToNextElm(qryTwoBitST);

        goto checkIfPrint;
      } // If have an insertion

      // Check what the reference postion is
      switch(getTwoBitAryElm(refTwoBitST))
      { // Switch; check if match, SNP, or deletion
        case defSnpFlag:
        // Case: Have a SNP
          *tmpRefCStr = *refSeqCStr;
          *tmpQryCStr = *qrySeqCStr;
          *tmpEqxCStr = 'X';

          ++refSeqCStr;
          ++qrySeqCStr;
          twoBitAryMoveToNextElm(refTwoBitST);
          twoBitAryMoveToNextElm(qryTwoBitST);
          break;
        // Case: Have a SNP

        case defMatchFlag:
        // Case: Have a Match
          *tmpRefCStr = *refSeqCStr;
          *tmpQryCStr = *qrySeqCStr;
          *tmpEqxCStr = '=';

          ++refSeqCStr;
          ++qrySeqCStr;
          twoBitAryMoveToNextElm(refTwoBitST);
          twoBitAryMoveToNextElm(qryTwoBitST);
          break;
        // Case: Have a Match
 
        case defGapFlag:
        // Case: have a deletion
          *tmpRefCStr = *refSeqCStr;
          *tmpQryCStr = '-';
          *tmpEqxCStr = 'D';

          ++refSeqCStr;
          twoBitAryMoveToNextElm(refTwoBitST);
          break;
        // Case: have a deletion
      } // Switch; check if match, SNP, or deletion

      checkIfPrint:

      if(lenBuffUL >= lnWrapUL - 1)
      { // IF at the end of the buffer
        // Write out the entry with the new line at the end
        fwrite(refAlnCStr,sizeof(char),lnWrapUL+1,alnFILE);
        fwrite(qryAlnCStr,sizeof(char),lnWrapUL+1,alnFILE);
        fwrite(eqxAlnCStr,sizeof(char),lnWrapUL+1,alnFILE);

        // Make sure a new line between entries
        fwrite("\n" , sizeof(char), 1, alnFILE);

        lenBuffUL = 0;
        tmpRefCStr = refAlnCStr;
        tmpQryCStr = qryAlnCStr;
        tmpEqxCStr = eqxAlnCStr;
      } // IF at the end of the buffer

      else
      { // Else I have more space in the buffer
        ++lenBuffUL;
        ++tmpRefCStr;
        ++tmpQryCStr;
        ++tmpEqxCStr;
      } // Else I have more space in the buffer
  } // While I have gaps at the start

  if(lenBuffUL > 0)
  { // If have a final bit of the buffer to print out
    // Write out the entry with the new line at the end
    fwrite(refAlnCStr, sizeof(char), lenBuffUL, alnFILE);
    fwrite("\n" , sizeof(char), 1, alnFILE);

    fwrite(qryAlnCStr, sizeof(char), lenBuffUL, alnFILE);
    fwrite("\n" , sizeof(char), 1, alnFILE);

    fwrite(eqxAlnCStr, sizeof(char), lenBuffUL, alnFILE);
    fwrite("\n" , sizeof(char), 1, alnFILE);
  } // If have a final bit of the buffer to print out

  free(refAlnCStr);
  refAlnCStr = 0;

  free(qryAlnCStr);
  qryAlnCStr = 0;

  free(eqxAlnCStr);
  eqxAlnCStr = 0;

  return 0;
} // printTwoBitAln
