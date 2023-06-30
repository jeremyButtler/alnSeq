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
  struct alnSet *settings, // Settings to use for alignment
  char useWaterScoreBl     // 1: Treat negatives as 0, like
                           // in a Smith Waterman
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

   struct twoBitAry *dirRowST = 0;  // Direction rows
   struct twoBitAry *twoBitAlnST = 0; // Holds alignment
   struct alnStruct *alignmentST = 0; // returned alignment

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

   twoBitAlnST = makeTwoBitArry(lenRefUL + lenQryUL, 0);
    // + 2 for index 0

   if(twoBitAlnST == 0) return 0;

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
     freeTwoBitAry(twoBitAlnST, 0, 0);
     return 0;
   } // If had a memory allocatoin error

   reverseScoreRowL = malloc(sizeof(long) * lenRefUL);

   if(reverseScoreRowL == 0)
   { // If had a memory allocatoin error
     freeTwoBitAry(twoBitAlnST, 0, 0);
     free(forwardScoreRowL);
     return 0;
   } // If had a memory allocatoin error

   /******************************************************\
   * Fun-01 Sec-02 Sub-03:
   *  - Initalize the direction rows
   \******************************************************/

   // I am using full length arrays to make the later
   // steps eaiser. This takes more memory, but makes life
   // nicer

   dirRowST = makeTwoBitArry(lenRefUL, 0);

   if(dirRowST == 0)
   { // If had a memory allocatoin error
     freeTwoBitAry(twoBitAlnST, 0, 0);
     free(forwardScoreRowL);
     free(reverseScoreRowL);
    
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
     dirRowST,        // For gap extension
     twoBitAlnST,     // Holds the alignment codes
     settings,        // Settings for the alignment
     useWaterScoreBl  // Tells if using negatives
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-04:
   ^    - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Mark the end of the alignment
   changeTwoBitElm(twoBitAlnST, defMoveStop);

   // Free the variables that are no longer used
   freeTwoBitAry(dirRowST, 0, 0);
   free(forwardScoreRowL);
   free(reverseScoreRowL);

   // Convert the alignment from a twoBitArray to an
   // Aligment structure. This step slows down the code,
   // but also reduces memory usage slightly
   alignmentST = twoBitAlnAryToAlnST(twoBitAlnST);
   freeTwoBitAry(twoBitAlnST, 0, 0);

   if(alignmentST == 0) return 0; // Memory error

   return alignmentST;
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
  unsigned long queryStartUL,//index 0 Starting query base
  unsigned long queryLenUL,//index 1 Length of target

  long *forwardScoreRowL,   // Holds final forward row
  long *reverseScoreRowL,   // For finding reverse scores
    // both the forward and reverse scoring rows must be
    // the size of the full length reference.
  struct twoBitAry *dirRowST,   // For gap extension
  struct twoBitAry *twoBitAlnST,// Holds alignment codes
  struct alnSet *settings, // Settings to use for alignment
  char useWaterScoreBl     // 1: Treat negatives as 0, like
                           // in a Smith Waterman
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
     for(
       unsigned long numInsUI = queryStartUL;
       numInsUI <= queryStartUL + queryLenUL;
       ++numInsUI
     ){ // Loop: fill in the insertions
       changeTwoBitElm(twoBitAlnST, defMoveLeft);
       twoBitAryMoveToNextElm(twoBitAlnST);
     } // Loop: fill in the insertions

     return; // Nothing else to do
   } // If all remaing bases are query insertions

   /******************************************************\
   * Fun-02 Sec-02 Sub-02:
   *  - Handle cases were I have just deletions
   \******************************************************/

   if(queryLenUL == 0)
   { // If all remaing bases are query deletions
     for(
       unsigned long numInsUI = refStartUL;
       numInsUI <= refStartUL + refLenUL;
       ++numInsUI
     ){ // Loop: fill in the insertions
       changeTwoBitElm(twoBitAlnST, defMoveUp);
       twoBitAryMoveToNextElm(twoBitAlnST);
     } // Loop: fill in the insertions

     return; // Nothing else to do
   } // If all remaing bases are query deletions

   /******************************************************\
   * Fun-02 Sec-02 Sub-03:
   *  - Handle cases were I have to align last ref base
   \******************************************************/

   if(refLenUL == 1)
   { // If I have to align the last reference base

     if(queryLenUL == 1)
     { // If bases are aligned (one reference & one query)
        changeTwoBitElm(twoBitAlnST, defMoveDiagnol);
        twoBitAryMoveToNextElm(twoBitAlnST);
        return; // Finished
     } // If bases are aligned (one reference & one query)

     positionSingleRefBase(
       *(refSeqCStr + refStartUL),// ref base
       qrySeqCStr + queryStartUL, // first base of query
       queryLenUL,                // Length of the query
       twoBitAlnST,              // Array to hold alignment
       settings,                  // Has Scoring variables
       useWaterScoreBl            // 0: keep negatives
     );

     return; // This base is now aligned
   } // If I have to align the last reference base

   /******************************************************\
   * Fun-02 Sec-02 Sub-04:
   *  - Handle cases were I have to align last query base
   \******************************************************/

   if(queryLenUL == 1)
   { // If I have to align the last query base

     if(refLenUL == 1)
     { // If bases are aligned (one query & one reference)
        changeTwoBitElm(twoBitAlnST, defMoveDiagnol);
        twoBitAryMoveToNextElm(twoBitAlnST);
        return; // Finished
     } // If bases are aligned (one query & one reference)

      // Get the scores, so that I can position the query
      // I do not care about the indel score
      scoreForwardHirsch(
        refSeqCStr,
        refStartUL,
        refLenUL,
        qrySeqCStr,         // Query seq with coordinates
        queryStartUL,
        queryLenUL,
        forwardScoreRowL,   // Array of scores to fill
        dirRowST,           // direction row for gap extend
        settings,           // setttings to use
        useWaterScoreBl     // 1: replace - numbers wihth 0
      );

     // Use the scoring array to fill in the alignment
     moveXElmFromStart(dirRowST, refStartUL + refLenUL);

     for(
       unsigned long endAlnUL = 0;
       endAlnUL >= refStartUL + refLenUL - 1;
       ++endAlnUL
     ) { // Loop: Fill in the direction matrix
       switch(getTwoBitAryElm(dirRowST))
       { // Check wich element I replace
         case defMoveLeft:
           changeTwoBitElm(twoBitAlnST, defMoveLeft);
           twoBitAryMoveBackOneElm(dirRowST);
           break;

         // At this point I found the base position
         // I also need to Make sure all future movements
         // are deletions
         case defMoveDiagnol: 
           changeTwoBitElm(twoBitAlnST, defMoveDiagnol);
           changeTwoBitElm(dirRowST, defMoveStop);
           break;

         case defMoveUp: 
           changeTwoBitElm(twoBitAlnST, defMoveUp);
           changeTwoBitElm(dirRowST, defMoveStop);
           break;

         // At this point I have flagged the remaining
         // bases as deletins
         case defMoveStop:
           changeTwoBitElm(twoBitAlnST, defMoveLeft);
       } // Check wich element I replace
         
       twoBitAryMoveToNextElm(twoBitAlnST);
     } // Loop: Fill in the direction matrix

     return; // Finshed aligning this query base
   } // If I have to align the last query base

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-03:
   ^  - Get scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    forwardIndelColL = 
      scoreForwardHirsch(
        refSeqCStr,       // Entire reference sequence
        refStartUL,       // Starting base of ref target
        refLenUL,         // length of ref target region
        qrySeqCStr,       // Query seq with coordinates
        queryStartUL,     // Starting base of query target
        queryLenUL / 2,   // Length of query target region
        forwardScoreRowL, // Array of scores to fill
        dirRowST,         // direction row for gap extend
        settings,         // setttings to use
        useWaterScoreBl   // 1: replace - numbers wihth 0
    ); // Get the scores for the forward direction

    reverseIndelColL = 
      scoreReverseHirsch(
        refSeqCStr,         // Entire reference sequence
        refStartUL,         // Starting base of ref target
        refLenUL,           // length of ref target region
        qrySeqCStr,         // Query seq with coordinates
        queryLenUL / 2,     // Start of query target
        queryLenUL / 2,     // Length  of query target
        forwardScoreRowL,   // Array of scores to fill
        dirRowST,           // direction row for gap extend
        settings,           // setttings to use
        useWaterScoreBl     // 1: replace - numbers wihth 0
      ); // Get the scores for the reverse direction
      // I can get away with queryLen/2 here, because 
      // queryLen is index 1 and the function takes in
      // an lenth 1 argument

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-04:
   ^   - Find the midpoint
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
     unsigned long lenRowUL = 0;
     lenRowUL < refLenUL + 2; // lenRefUI is index 0
     ++lenRowUL
   ) { // Loop; add up all scores
     *(forwardScoreRowL + lenRowUL) += 
       *(reverseScoreRowL + lenRowUL);
       // The reverse row is already reversed

     // Check if ths midpoint is better
     if(
       *(forwardScoreRowL + lenRowUL)
         > *(forwardScoreRowL + midPointUL)
     ) midPointUL = lenRowUL;
   }// Loop; add up all scores

   forwardIndelColL += reverseIndelColL;

   if(*(forwardScoreRowL + midPointUL) < forwardIndelColL)
     midPointUL = 0;
   else ++midPointUL; // Convert to length (index 1)

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-05:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   HirschbergFun(
     refSeqCStr,      // Full reference sequence
     refStartUL,      // Full reference sequence
     midPointUL,
     qrySeqCStr,      // Full queyr sequence
     queryStartUL,    // Start of queyr target region
     queryLenUL / 2,  // Length of query target region
     forwardScoreRowL,  // For scoring
     reverseScoreRowL, // Has last line of scores
     dirRowST,           // direction row for gap extend
     twoBitAlnST,     // Holds the alignment codes
     settings,        // Settings for the alignment
     useWaterScoreBl  // Tells if using negatives
   );

   HirschbergFun(
     refSeqCStr,                 // Full reference sequence
     refStartUL + midPointUL - 1,  // New reference start
     refLenUL - midPointUL,        // New reference end
     qrySeqCStr,                 // Full query sequence
     queryStartUL + (queryLenUL / 2) + 1,
     (queryLenUL / 2) - 1,       // New query length
     forwardScoreRowL,  // For scoring
     reverseScoreRowL, // Has last line of scores
     dirRowST,           // direction row for gap extend
     twoBitAlnST,     // Holds the alignment codes
     settings,        // Settings for the alignment
     useWaterScoreBl  // Tells if using negatives
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
  unsigned long queryStartUL,// Index 0 Starting query base
  unsigned long queryLenUL,  // index 1 length of target

  long *scoreRowPtrL,        // Array of scores to fill
  struct twoBitAry *dirRowST,//direction row for gap extend
  struct alnSet *settings,   // setttings to use
  char useWaterScoreBl       //1: replace - numbers wihth 0
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

   long *scoreOnPtrL = scoreRowPtrL;
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

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1: break; // Already 0

     case 0:
     default:
       *scoreOnPtrL = settings->gapStartPenaltyI;
       changeTwoBitElm(dirRowST, defMoveLeft);
       break;
   } // Switch: Check if changing negatives to 0's

   ++scoreOnPtrL;
   twoBitAryMoveToNextElm(dirRowST);

   /******************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Set up remaing elements (indels) (uses gap extend)
   \******************************************************/

   for(
     char *refBaseCStr = refSeqCStr + refStartUL;
     refBaseCStr < refSeqCStr + refStartUL + refLenUL;
     ++refBaseCStr
   ){ // Loop:Set the initial blank scores
       switch(useWaterScoreBl)
       { // Switch: Check if changing negatives to 0's
         case 1:
           *scoreOnPtrL = 0;
           changeTwoBitElm(dirRowST, defMoveStop);
           break;

         case 0:
         default:
           *scoreOnPtrL =
             *(scoreOnPtrL-1) +settings->gapExtendPenaltyI;

           changeTwoBitElm(dirRowST, defMoveLeft);
           break;
       } // Switch: Check if changing negatives to 0's

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
   scoreOnPtrL = scoreRowPtrL;

   // Fine the first match/snp score (first ref base)
   nextMatchScoreL =
     getBasePairScore(
       qrySeqCStr + queryStartUL,   // first query base
       refSeqCStr + refStartUL,     // first ref base
       settings                     // Has score matrix
   ); // Ge the score for the enxt base

   nextMatchScoreL += *scoreOnPtrL;
     // This is here so I can overwrite the array with the
     // new scores

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1: break; // indelColL already 0

     case 0:
     default:
       indelColL += settings->gapStartPenaltyI;
       break;
   } // Switch: Check if changing negatives to 0's

   // Find the first insertion and deletion scores
   insScoreL = indelColL + settings->gapExtendPenaltyI;
   delScoreL = indelColL + settings->gapExtendPenaltyI;

   /******************************************************\
   * Fun-03 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   for(
     char *qryBaseCStr = qrySeqCStr + queryStartUL;
     qryBaseCStr < qrySeqCStr + queryStartUL + queryLenUL;
     ++qryBaseCStr
   ){ // Loop: score all query bases (rows)

     for(
       char *refBaseCStr = refSeqCStr + refStartUL;
       refBaseCStr < refSeqCStr + refStartUL + refLenUL;
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

       // Set your score and direction
       switch(useWaterScoreBl)
       { // Switch: check if making 0's
         case 1:
           updateDirAndScoreNeedle(
             dirRowST,        // Direction matrix
             settings,        // has direction preference
             &insScoreL,      // Score for an insertion
             &matchScoreL,    // Score for an deletion
             &delScoreL,      // Score for an match/snp
             scoreOnPtrL        // Score position to update
           ); // Update the score using a waterman
           break;

         case 0:
         default:
           updateDirAndScoreNeedle(
             dirRowST,        // Direction matrix
             settings,        // has direction preference
             &insScoreL,      // Score for an insertion
             &matchScoreL,    // Score for an deletion
             &delScoreL,      // Score for an match/snp
             scoreOnPtrL        // Score position to update
           ); // Update the scores using a needleman

           break;
       } // Switch: check if making 0's

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

     // Move back to the start of row
     moveXElmFromStart(dirRowST, refStartUL);
     scoreOnPtrL = scoreRowPtrL;

     // Fine the first match/snp score (first ref base)
     nextMatchScoreL =
       getBasePairScore(
         qryBaseCStr + 1,           // Next query base
         refSeqCStr + refStartUL,   // first base
         settings                   // Has score matrix
     ); // Ge the score for the enxt base

     nextMatchScoreL += *scoreOnPtrL;
       // This is here so I can overwrite the array with the
       // new scores

     switch(useWaterScoreBl)
     { // Switch: Check if changing netatives 0's
       case 1: break; // indelColL already 0

       case 0:
       default:
         indelColL += settings->gapStartPenaltyI;
         break;
     } // Switch: Check if changing negatives to 0's
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
  unsigned long queryStartUL,// Index 0 Starting query base
  unsigned long queryLenUL,  // index 1 length of target

  long *scoreRowPtrL,       // Array of scores to fill
     // This needs to be as long as the full length
     // reference sequence
  struct twoBitAry *dirRowST,//direction row for gap extend
  struct alnSet *settings,  // setttings to use
  char useWaterScoreBl      //1: replace - numbers wihth 0
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

   long *scoreOnPtrL = scoreRowPtrL + refLenUL - 1;
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

   moveXElmFromStart(dirRowST, refLenUL - 1);
     // - 1 to account for refLenUL being index 1


   indelColL = 0;

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1: break; // already 0

     case 0:
     default:
       *scoreOnPtrL = settings->gapStartPenaltyI;
       changeTwoBitElm(dirRowST, defMoveLeft);
       break;
   } // Switch: Check if changing negatives to 0's

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
       refSeqCStr + refStartUL + refLenUL - 1;
     refBaseCStr > refSeqCStr + refStartUL - 1; // index 0
     --refBaseCStr
   ){ // Loop:Set the initial blank scores
       switch(useWaterScoreBl)
       { // Switch: Check if changing negatives to 0's
         case 1:
           *scoreOnPtrL = 0;
           changeTwoBitElm(dirRowST, defMoveStop);
           break;

         case 0:
         default:
           *scoreOnPtrL =
             *(scoreOnPtrL+1) +settings->gapExtendPenaltyI;

           changeTwoBitElm(dirRowST, defMoveLeft);
           break;
       } // Switch: Check if changing negatives to 0's

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
   moveXElmFromStart(dirRowST, refLenUL - 1);
   scoreOnPtrL = scoreRowPtrL + refLenUL - 1;

   // Fine the first match/snp score (first ref base)
   nextMatchScoreL =
     getBasePairScore(
       qrySeqCStr + queryStartUL + queryLenUL - 1,
       refSeqCStr + refStartUL + refLenUL - 1,
       settings                         // Has score matrix
   ); // Ge the score for the enxt base

   nextMatchScoreL += *scoreOnPtrL;
     // This is here so I can overwrite the array with the
     // new scores

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1: break; // indelColL already 0

     case 0:
     default:
       indelColL += settings->gapStartPenaltyI;
       break;
   } // Switch: Check if changing negatives to 0's

   // Find the first insertion and deletion scores
   insScoreL = indelColL + settings->gapExtendPenaltyI;
   delScoreL = indelColL + settings->gapExtendPenaltyI;

   /******************************************************\
   * Fun-04 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   for(
     char *qryBaseCStr =
       qrySeqCStr + queryStartUL + queryLenUL - 1;
     qryBaseCStr > qrySeqCStr + queryStartUL - 1;
     --qryBaseCStr
   ){ // Loop: score all query bases (rows)

     for(
       char *refBaseCStr =
         refSeqCStr + refStartUL + refLenUL - 1;
       refBaseCStr > refSeqCStr + refStartUL - 1;
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

       // Set your score and direction
       switch(useWaterScoreBl)
       { // Switch: check if making 0's
         case 1:
           updateDirAndScoreNeedle(
             dirRowST,        // Direction matrix
             settings,        // has direction preference
             &insScoreL,      // Score for an insertion
             &matchScoreL,    // Score for an deletion
             &delScoreL,      // Score for an match/snp
             scoreOnPtrL        // Score position to update
           ); // Update the score using a waterman
           break;

         case 0:
         default:
           updateDirAndScoreNeedle(
             dirRowST,        // Direction matrix
             settings,        // has direction preference
             &insScoreL,      // Score for an insertion
             &matchScoreL,    // Score for an deletion
             &delScoreL,      // Score for an match/snp
             scoreOnPtrL        // Score position to update
           ); // Update the scores using a needleman

           break;
       } // Switch: check if making 0's

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

     // Move back to the start of row
     moveXElmFromStart(dirRowST, refStartUL + refLenUL -1);
     scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL -1;

     // Fine the first match/snp score (first ref base)
     nextMatchScoreL =
       getBasePairScore(
         qryBaseCStr - 1,           // Next query base
         refSeqCStr + refStartUL + refLenUL - 1,
         settings         // Has score matrix
     ); // Ge the score for the enxt base

     nextMatchScoreL += *scoreOnPtrL;
       // This is here so I can overwrite the array with the
       // new scores

     switch(useWaterScoreBl)
     { // Switch: Check if changing netatives 0's
       case 1: break; // indelColL already 0

       case 0:
       default:
         indelColL += settings->gapStartPenaltyI;
         break;
     } // Switch: Check if changing negatives to 0's

     --scoreOnPtrL;
     twoBitAryMoveBackOneElm(dirRowST);
     // the lastElm is currently at the correct position
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
|    o twoBitAlnST to hold the aligned reference and query
|      bases
\--------------------------------------------------------*/
void positionSingleRefBase(
  char refBaseC,          // Single reference base to align
  char *qrySeqCStr,//query sequence to position ref base on
  unsigned long endOfQrySeqUL,   // Marks end of alignment
  struct twoBitAry *twoBitAlnST, // Array to hold alignment
  struct alnSet *settings,       // setttings to use
  char useWaterScoreBl          // convert negatives to 0's
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBase
   '  - Align a single reference base to a query sequence
   '    This is a really odd situation, were I have only
   '    the indel column and  single base in a row, but
   '    also have a large query.
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
   long indelColScoreL = 0; // Current score to add to del
   long lastIndelCoL = 0;   // Score to add to a match

   char *qryBaseCStr = qrySeqCStr;
   char *lastMatchCStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-02:
   ^  - Find the reference bases position on the query
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   changeTwoBitElm(twoBitAlnST, defMoveStop);
   indelColScoreL += settings->gapStartPenaltyI;

   while(qryBaseCStr <= qrySeqCStr + endOfQrySeqUL)
   { // While I have bases to compare
     insScoreL =
       getIndelScore(twoBitAlnST, settings, &curScoreL);

     delScoreL =
       indelColScoreL + settings->gapExtendPenaltyI;

     matchScoreL = lastIndelCoL;

     matchScoreL +=
       getBasePairScore(qryBaseCStr, &refBaseC, settings);

     // Set your score and direction
     switch(useWaterScoreBl)
     { // Switch: check if making 0's
       case 1:
         updateDirAndScoreNeedle(
           twoBitAlnST,     // Direction matrix
           settings,        // has direction preference
           &insScoreL,      // Score for an insertion
           &matchScoreL,    // Score for an deletion
           &delScoreL,      // Score for an match/snp
           &curScoreL        // Score position to update
         ); // Update the score using a waterman
         break;

       case 0:
       default:
         updateDirAndScoreNeedle(
           twoBitAlnST,     // Direction matrix
           settings,        // has direction preference
           &insScoreL,      // Score for an insertion
           &matchScoreL,    // Score for an deletion
           &delScoreL,      // Score for an match/snp
           &curScoreL        // Score position to update
         ); // Update the scores using a needleman

         break;
     } // Switch: check if making 0's

     switch(getTwoBitAryElm(twoBitAlnST))
     { // Switch: Check if keeping the score
       case defMoveDiagnol:
         lastMatchCStr = qryBaseCStr;
         break;

       case defMoveLeft: break;
       case defMoveUp: break;
       case defMoveStop: break;
     } // Switch: Check if keeping the score

     // Update indel column score and find deletion score
     lastIndelCoL = indelColScoreL;
     indelColScoreL += settings->gapExtendPenaltyI;

     ++qryBaseCStr;
   } // While I have bases to compare

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-03:
   ^  - Fill in the insertions and reference base position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   qryBaseCStr = qrySeqCStr;

   // Add in the insertions at the start
   while(qryBaseCStr < lastMatchCStr)
   { // While I have insertions to fill
     changeTwoBitElm(twoBitAlnST, defMoveUp);
     twoBitAryMoveToNextElm(twoBitAlnST);
     ++qryBaseCStr;
   } // While I have insertions to fill
   
   // Add in the positoin of the base
   changeTwoBitElm(twoBitAlnST, defMoveDiagnol);
   twoBitAryMoveToNextElm(twoBitAlnST);

   // Finish adding in the insertions at the end
   while(qryBaseCStr <= qrySeqCStr + endOfQrySeqUL)
   { // While I have insertions to fill
     changeTwoBitElm(twoBitAlnST, defMoveUp);
     twoBitAryMoveToNextElm(twoBitAlnST);
     ++qryBaseCStr;
   } // While I have insertions to fill

   return;
} // positionSingleRefBase

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o An alnStruct with the alignment from twoBitAlnST
|  - Modifies:
|    o twoBitAlnST to point to the end of the alignment
\--------------------------------------------------------*/
struct alnStruct * twoBitAlnAryToAlnST(
  struct twoBitAry *twoBitAlnST
   // Two bit array with alignment to convert, which ends
   // in a stop
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: twoBitAlnAryToAlnST
   '  - Converts a two bit array with an alignment to an
   '    alnStruct structure for printing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   uint8_t dirElmUC = 0; // Direction to convert to
   struct alnStruct *alnST = 0;

   // Set up the alignment strucuter
   alnST = malloc(sizeof(struct alnStruct));

   if(alnST == 0) return 0;

   initAlnST(alnST);

  if(addAlnSTArray(alnST, lenTwoBitAry(twoBitAlnST)) != 0)
  { // If had a memory allocation error
    freeAlnST(alnST, 1);
    return 0;
  } // If had a memory allocation error

  // Move to the start of the array and get first element
  moveXElmFromStart(twoBitAlnST, 0);
  dirElmUC = getTwoBitAryElm(twoBitAlnST);

  // Add the twobit array directions to alignment structure
  while(dirElmUC != defMoveStop)
  { // While I have array elements to convert
    alnSTAddNewCode(alnST, dirElmUC);
    twoBitAryMoveToNextElm(twoBitAlnST);
    dirElmUC = getTwoBitAryElm(twoBitAlnST);
  } // While I have array elements to convert

  return alnST;
} // twoBitAlnAryToAlnST
