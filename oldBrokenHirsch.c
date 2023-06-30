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
    // For refST and qryST, use seqStruct->offsetUI to set
    // the starting point for the alignmnet and
    // seqStruct->endAlnUI to set the ending point
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
   
   uint32_t lenRefUI = refST->endAlnUI - refST->offsetUI;
   uint32_t lenQryUI = qryST->endAlnUI - qryST->offsetUI;

   uint32_t refOffsetUI = refST->offsetUI;
   uint32_t refEndAlnUI = refST->endAlnUI;

   uint32_t qryOffsetUI = qryST->offsetUI;
   uint32_t qryEndAlnUI = qryST->endAlnUI;

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

   twoBitAlnST = makeTwoBitArry(lenRefUI + lenQryUI +2, 0);
    // + 2 for index 0

   if(twoBitAlnST == 0) return 0;

   /******************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Initalize the scoring rows
   \******************************************************/

   // I am using full length arrays to make the later
   // steps eaiser. This takes more memory, but makes life
   // nicer

   forwardScoreRowL = malloc(sizeof(long) * (lenRefUI+2));

   if(forwardScoreRowL == 0)
   { // If had a memory allocatoin error
     freeTwoBitAry(twoBitAlnST, 0, 0);
     return 0;
   } // If had a memory allocatoin error

   reverseScoreRowL = malloc(sizeof(long) * (lenRefUI+2));

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

   dirRowST = makeTwoBitArry(lenRefUI + 2, 0);

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
     refST,           // Has the start and end coordianates
     qryST,
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

   // Rest the alignment values
   refST->offsetUI = refOffsetUI;
   refST->endAlnUI = refEndAlnUI;

   qryST->offsetUI = qryOffsetUI;
   qryST->endAlnUI = qryEndAlnUI;

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
  struct seqStruct *refST, // Reference sequence to align
  struct seqStruct *qryST, // Qeury sequence to align
    // For refST and qryST, use seqStruct->offsetUI to set
    // the starting point for the alignmnet and
    // seqStruct->endAlnUI to set the ending point
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

   uint32_t lenRefUI = refST->endAlnUI - refST->offsetUI;
   uint32_t lenQryUI = qryST->endAlnUI - qryST->offsetUI;

   uint32_t refOffsetUI = refST->offsetUI;
   uint32_t refEndAlnUI = refST->endAlnUI;

   uint32_t qryOffsetUI = qryST->offsetUI;
   uint32_t qryEndAlnUI = qryST->endAlnUI;

   uint32_t midPointUI = 0;

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

   if(lenRefUI == 0)
   { // If all remaing bases are query insertions
     for(
       uint32_t numInsUI = qryST->offsetUI;
       numInsUI <= qryST->endAlnUI;
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

   if(lenQryUI == 0)
   { // If all remaing bases are query deletions
     for(
       uint32_t numInsUI = refST->offsetUI;
       numInsUI <= refST->endAlnUI;
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

   if(lenRefUI == 1)
   { // If I have to align the last reference base

     if(lenQryUI == 1)
     { // If bases are aligned (one reference & one query)
        changeTwoBitElm(twoBitAlnST, defMoveDiagnol);
        twoBitAryMoveToNextElm(twoBitAlnST);
        return; // Finished
     } // If bases are aligned (one reference & one query)

     positionSingleRefBase(
       *(refST->seqCStr + refST->offsetUI),
       qryST->seqCStr + qryST->offsetUI,
       qryST->endAlnUI - qryST->offsetUI,
       twoBitAlnST, // Array to hold alignment
       settings,     // setttings to use
       useWaterScoreBl
     );

     return; // This base is now aligned
   } // If I have to align the last reference base

   /******************************************************\
   * Fun-02 Sec-02 Sub-04:
   *  - Handle cases were I have to align last query base
   \******************************************************/

   if(lenQryUI == 1)
   { // If I have to align the last query base

     if(lenRefUI == 1)
     { // If bases are aligned (one query & one reference)
        changeTwoBitElm(twoBitAlnST, defMoveDiagnol);
        twoBitAryMoveToNextElm(twoBitAlnST);
        return; // Finished
     } // If bases are aligned (one query & one reference)

      // Get the scores, so that I can position the query
      scoreForwardHirsch(
        refST,
        qryST,              // Query seq with coordinates
        forwardScoreRowL,   // Array of scores to fill
        dirRowST,           // direction row for gap extend
        settings,           // setttings to use
        useWaterScoreBl     // 1: replace - numbers wihth 0
      );

     // Use the scoring array to fill in the alignment
     moveXElmFromStart(dirRowST, refST->endAlnUI);

     for(
       uint32_t endAlnUI = refST->endAlnUI;
       endAlnUI >= refST->endAlnUI;
       --endAlnUI
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

    // Else
    qryST->endAlnUI = qryOffsetUI + (qryEndAlnUI / 2);
    qryST->offsetUI = qryOffsetUI;

    scoreForwardHirsch(
      refST,
      qryST,              // Query seq with coordinates
      forwardScoreRowL,   // Array of scores to fill
      dirRowST,           // direction row for gap extend
      settings,           // setttings to use
      useWaterScoreBl     // 1: replace - numbers wihth 0
    );

    qryST->offsetUI = qryST->endAlnUI + 1;
    qryST->endAlnUI = qryEndAlnUI;

    scoreReverseHirsch(
      refST,
      qryST,              // Query seq with coordinates
      reverseScoreRowL,   // Array of scores to fill
      dirRowST,           // direction row for gap extend
      settings,           // setttings to use
      useWaterScoreBl     // 1: replace - numbers wihth 0
    );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-04:
   ^   - Find the midpoint
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
     uint32_t lenRowUI = 0;
     lenRowUI < lenRefUI + 2; // lenRefUI is index 0
     ++lenRowUI
   ) { // Loop; add up all scores
     *(forwardScoreRowL + lenRowUI) += 
       *(reverseScoreRowL + lenRowUI);
       // The reverse row is already reversed

     // Check if ths midpoint is better
     if(
       *(forwardScoreRowL + lenRowUI)
         > *(forwardScoreRowL + midPointUI)
     ) midPointUI = lenRowUI;
   }// Loop; add up all scores

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-05:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   qryST->offsetUI = qryOffsetUI;
   qryST->endAlnUI = qryOffsetUI + (lenQryUI / 2);
     // lenQryUI is index 1

   // This is not the best way, but it works
   // Handle case of all bases go to next query half
   // This is because I am working with 0 index
   if(midPointUI == 0) refST->endAlnUI = 0;
   else refST->endAlnUI = midPointUI - 1 + refOffsetUI;

   refST->offsetUI = refOffsetUI;
       // -1 to handle the indel column

   HirschbergFun(
     refST,           // Has the start and end coordianates
     qryST,
     forwardScoreRowL,  // For scoring
     reverseScoreRowL, // Has last line of scores
     dirRowST,           // direction row for gap extend
     twoBitAlnST,     // Holds the alignment codes
     settings,        // Settings for the alignment
     useWaterScoreBl  // Tells if using negatives
   );

   qryST->offsetUI = qryOffsetUI + 1 + (lenQryUI /2);
   qryST->endAlnUI = qryEndAlnUI;

   // Handle case of all bases go to previous query half
   // This is because I am working with 0 index
   if(midPointUI > lenRefUI) refST->offsetUI = refEndAlnUI;
   else if(midPointUI == lenRefUI)
     refST->offsetUI = refEndAlnUI - 1;
   else refST->offsetUI = midPointUI + refOffsetUI;// + 1;

   refST->endAlnUI = refEndAlnUI;

   HirschbergFun(
     refST,           // Has the start and end coordianates
     qryST,
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

   // Rest the alignment values
   refST->offsetUI = refOffsetUI;
   refST->endAlnUI = refEndAlnUI;

   qryST->offsetUI = qryOffsetUI;
   qryST->endAlnUI = qryEndAlnUI;

   return;
} // HirschbergFun

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      Needleman Wunsch / Smith Waterman alignment
|    o dirRowSt to hold the last row of directions in a
|      Needleman Wunsch / Smith Waterman alignment
\--------------------------------------------------------*/
void scoreForwardHirsch(
  struct seqStruct *refST,   // Ref seq with coordinates
  struct seqStruct *qryST,   // Query seq with coordinates
    // Coodanitates:
    //  Starting: seqStruct->offsetUI (index 0)
    //  Ending: seqStruct->endAlnUI   (index 0)
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
   long insScoreL = 0;
   long delScoreL = 0;
   long matchScoreL = 0;
   long nextMatchScoreL = 0;

   struct twoBitAry lastElm;

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

   moveXElmFromStart(dirRowST, 0);

   *scoreOnPtrL = 0;
   changeTwoBitElm(dirRowST, defMoveStop);

   ++scoreOnPtrL;
   twoBitAryMoveToNextElm(dirRowST);

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1:
       *scoreOnPtrL = 0;
       changeTwoBitElm(dirRowST, defMoveStop);
       break;

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
     char *refBaseCStr = refST->seqCStr + refST->offsetUI;
     refBaseCStr < refST->seqCStr + refST->endAlnUI + 1;
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
   moveXElmFromStart(dirRowST, 0);
   cpTwoBitPos(dirRowST, &lastElm);
   scoreOnPtrL = scoreRowPtrL;

   // Fine the first match/snp score (first ref base)
   nextMatchScoreL =
     getBasePairScore(
       qryST->seqCStr + qryST->offsetUI,// first query base
       refST->seqCStr + refST->offsetUI,// first ref base
       settings                         // Has score matrix
   ); // Ge the score for the enxt base

   nextMatchScoreL += *scoreOnPtrL;
     // This is here so I can overwrite the array with the
     // new scores

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1:
       *scoreOnPtrL = 0;
       changeTwoBitElm(dirRowST, defMoveStop);
       break;

     case 0:
     default:
       *scoreOnPtrL =
         *scoreOnPtrL + settings->gapStartPenaltyI;

       changeTwoBitElm(dirRowST, defMoveUp);
       break;
   } // Switch: Check if changing negatives to 0's

   ++scoreOnPtrL;
   twoBitAryMoveToNextElm(dirRowST);
   // the lastElm is currently at the correct position

   /******************************************************\
   * Fun-03 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   for(
     char *qryBaseCStr = qryST->seqCStr + qryST->offsetUI;
     qryBaseCStr < qryST->seqCStr + qryST->endAlnUI + 1;
     ++qryBaseCStr
   ){ // Loop: score all query bases (rows)

     for(
       char *refBaseCStr = refST->seqCStr +refST->offsetUI;
       refBaseCStr < refST->seqCStr + refST->endAlnUI + 1;
       ++refBaseCStr
     ){ // Loop:score all query bases (columns)

       // Get the insertion score
       insScoreL =
         getIndelScore(dirRowST, settings, scoreOnPtrL);

       // Get the deletion score
       delScoreL =
         getIndelScore(&lastElm, settings, scoreOnPtrL -1);
         // This is fun-03 in generalAlnFun.h

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

       twoBitAryMoveToNextElm(dirRowST);
       twoBitAryMoveToNextElm(&lastElm);
       ++scoreOnPtrL;
     } // Loop:score all query bases (columns)

     /****************************************************\
     * Fun-03 Sec-03 Sub-04:
     *  - Set up for scoring next row (handle indel col)
     \****************************************************/

     // Move back to the start of row
     moveXElmFromStart(dirRowST, 0);
     moveXElmFromStart(&lastElm, 0);
     scoreOnPtrL = scoreRowPtrL;

     // Fine the first match/snp score (first ref base)
     nextMatchScoreL =
       getBasePairScore(
         qryBaseCStr + 1,           // Next query base
         refST->seqCStr + refST->offsetUI, // first base
         settings         // Has score matrix
     ); // Ge the score for the enxt base

     nextMatchScoreL += *scoreOnPtrL;
       // This is here so I can overwrite the array with the
       // new scores

     switch(useWaterScoreBl)
     { // Switch: Check if changing netatives 0's
       case 1:                 // waterman setting (no -'s)
         *scoreOnPtrL = 0;
         changeTwoBitElm(dirRowST, defMoveStop);
         break;

       case 0:                 // neeldeman setting
       default:
         *scoreOnPtrL =
           *scoreOnPtrL + settings->gapExtendPenaltyI;

         changeTwoBitElm(dirRowST, defMoveUp);
         break;
     } // Switch: Check if changing negatives to 0's

     ++scoreOnPtrL;
     twoBitAryMoveToNextElm(dirRowST);
     // the lastElm is currently at the correct position

   } // Loop: score all query bases (rows)
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Clean up
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return;
} // scoreForwardHirsch

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      backwards Needleman Wunsch /Smith Waterman alignment
|    o dirRowST to hold the last row of directions in a
|      backwards Needleman Wunsch /Smith Waterman alignment
\--------------------------------------------------------*/
void scoreReverseHirsch(
  struct seqStruct *refST,   // Ref seq with coordinates
  struct seqStruct *qryST,   // Query seq with coordinates
    // Coodanitates:
    //  Starting: seqStruct->offsetUI (index 0)
    //  Ending: seqStruct->endAlnUI   (index 0)
  long *scoreRowPtrL,        // Array of scores to fill
     // This needs to be as long as the full length
     // reference sequence
  struct twoBitAry *dirRowST,//direction row for gap extend
  struct alnSet *settings,   // setttings to use
  char useWaterScoreBl       //1: replace - numbers wihth 0
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

   long lenRefUI = refST->endAlnUI - refST->offsetUI;
   long *scoreOnPtrL = scoreRowPtrL + lenRefUI + 1;
     // +1 to account for the indel column

   long insScoreL = 0;
   long delScoreL = 0;
   long matchScoreL = 0;
   long nextMatchScoreL = 0;

   struct twoBitAry lastElm;

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

   moveXElmFromStart(dirRowST, lenRefUI + 1);
     // +1 to account for the indel column

   *scoreOnPtrL = 0;
   changeTwoBitElm(dirRowST, defMoveStop);

   --scoreOnPtrL;
   twoBitAryMoveBackOneElm(dirRowST);

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1:
       *scoreOnPtrL = 0;
       changeTwoBitElm(dirRowST, defMoveStop);
       break;

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
     char *refBaseCStr = refST->seqCStr +refST->endAlnUI-1;
     refBaseCStr > refST->seqCStr + refST->offsetUI - 1;
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

   // Move back to the start of row
   moveXElmFromStart(dirRowST, lenRefUI + 1);
   cpTwoBitPos(dirRowST, &lastElm);
   scoreOnPtrL = scoreRowPtrL + lenRefUI + 1;

   // Fine the first match/snp score (first ref base)
   nextMatchScoreL =
     getBasePairScore(
       qryST->seqCStr + qryST->endAlnUI,
       refST->seqCStr + refST->endAlnUI,
       settings                         // Has score matrix
   ); // Ge the score for the enxt base

   nextMatchScoreL += *scoreOnPtrL;
     // This is here so I can overwrite the array with the
     // new scores

   switch(useWaterScoreBl)
   { // Switch: Check if changing negatives to 0's
     case 1:
       *scoreOnPtrL = 0;
       changeTwoBitElm(dirRowST, defMoveStop);
       break;

     case 0:
     default:
       *scoreOnPtrL =
         *scoreOnPtrL + settings->gapStartPenaltyI;

       changeTwoBitElm(dirRowST, defMoveUp);
       break;
   } // Switch: Check if changing negatives to 0's

   --scoreOnPtrL;
   twoBitAryMoveBackOneElm(dirRowST);
   // the lastElm is currently at the correct position

   /******************************************************\
   * Fun-04 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   for(
     char *qryBaseCStr = qryST->seqCStr + qryST->endAlnUI;
     qryBaseCStr > qryST->seqCStr + qryST->offsetUI - 1;
     --qryBaseCStr
   ){ // Loop: score all query bases (rows)

     for(
       char *refBaseCStr= refST->seqCStr + refST->endAlnUI;
       refBaseCStr > refST->seqCStr + refST->offsetUI - 1;
       --refBaseCStr
     ){ // Loop:score all query bases (columns)

       // Get the insertion score
       insScoreL =
         getIndelScore(dirRowST, settings, scoreOnPtrL);

       // Get the deletion score
       delScoreL =
         getIndelScore(&lastElm, settings, scoreOnPtrL +1);
         // This is fun-03 in generalAlnFun.h

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

       twoBitAryMoveBackOneElm(dirRowST);
       twoBitAryMoveBackOneElm(&lastElm);
       --scoreOnPtrL;
     } // Loop:score all query bases (columns)

     /****************************************************\
     * Fun-04 Sec-03 Sub-04:
     *  - Set up for scoring next row (handle indel col)
     \****************************************************/

     // Move back to the start of row
     moveXElmFromStart(dirRowST, lenRefUI + 1);
     moveXElmFromStart(&lastElm, lenRefUI + 1);
     scoreOnPtrL = scoreRowPtrL + lenRefUI + 1;

     // Fine the first match/snp score (first ref base)
     nextMatchScoreL =
       getBasePairScore(
         qryBaseCStr - 1,           // Next query base
         refST->seqCStr + refST->endAlnUI,
         settings         // Has score matrix
     ); // Ge the score for the enxt base

     nextMatchScoreL += *scoreOnPtrL;
       // This is here so I can overwrite the array with the
       // new scores

     switch(useWaterScoreBl)
     { // Switch: Check if changing netatives 0's
       case 1:                 // waterman setting (no -'s)
         *scoreOnPtrL = 0;
         changeTwoBitElm(dirRowST, defMoveStop);
         break;

       case 0:                 // neeldeman setting
       default:
         *scoreOnPtrL =
           *scoreOnPtrL + settings->gapExtendPenaltyI;

         changeTwoBitElm(dirRowST, defMoveUp);
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

   return;
} // scoreReverseHirsch

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the aligned reference and query
|      bases
\--------------------------------------------------------*/
void positionSingleRefBase(
  char refBaseC,    // Single reference base to align
  char *qrySeqCStr,//query sequence to position ref base on
  uint32_t endOfQrySeqUI,       // Marks end of alignment
  struct twoBitAry *twoBitAlnST, // Array to hold alignment
  struct alnSet *settings,      // setttings to use
  char useWaterScoreBl  // convert negatives to 0's
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

   while(qryBaseCStr <= qrySeqCStr + endOfQrySeqUI)
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
   while(qryBaseCStr <= qrySeqCStr + endOfQrySeqUI)
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
