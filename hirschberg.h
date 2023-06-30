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

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
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

#ifndef HIRSCHBERG_H
#define HIRSCHBERG_H

#include "needleman.h"
#include "waterman.h"

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o A twoBitAry structure with the alignment.
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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

  long *scoreRowPtrL,        // Array of scores to fill
     // This needs to be as long as the full length
     // reference sequence
  struct twoBitAry *dirRowST,//direction row for gap extend
  struct alnSet *settings,   // setttings to use
  char useWaterScoreBl       //1: replace - numbers wihth 0
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the aligned reference and query
|      bases
\--------------------------------------------------------*/
void positionSingleRefBase(
  char refBaseC,          // Single reference base to align
  char *qrySeqCStr,//query sequence to position ref base on
  unsigned long endOfQrySeqUI,   // Marks end of alignment
  struct twoBitAry *twoBitAlnST, // Array to hold alignment
  struct alnSet *settings,       // setttings to use
  char useWaterScoreBl          // convert negatives to 0's
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBase
   '  - Align a single reference base to a query sequence
   '  o fun-05 sec-01:
   '    - Variable declerations
   '  o fun-05 sec-02:
   '    - Find the reference bases position on the query
   '  o fun-05 sec-03:
   '    - Fill in insertions and reference base position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: twoBitAlnAryToAlnST
   '  - Converts a two bit array with an alignment to an
   '    alnStruct structure for printing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
