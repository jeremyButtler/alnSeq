/*#########################################################
# Name: hirschberg
# Use:
#  - Holds functions for doing a hirschberg global
#    alignment
# Libraries:
#  - needleman.h
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
  struct alnSet *settings  // Settings to use for alignment
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

  struct twoBitAry *refAlnST,   // For gap extension
  struct twoBitAry *qryAlnST,  // Holds alignment codes
  struct alnSet *settings  // Settings to use for alignment
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
  struct alnSet *settings    // setttings to use
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
  struct alnSet *settings    // setttings to use
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBase
   '  - Align a single base to a sequence
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: printTwoBitAln
   '  - Converts a two bit array with an alignment to an
   '    alnStruct structure for printing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
