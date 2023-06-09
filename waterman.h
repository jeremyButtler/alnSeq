/*########################################################
# Name waterman
# Use:
#  o Holds functions doing a Waterman-Smith pairwise
#    alignments. This version outputs a single alignment
# Includes:
#   - "generalAlnFun.h"
#   - "alnStruct.h"
#   - "alnMatrixStruct.h"
#   o "twoBitArrays.h"
#   o "scoresST.h"
#   o "seqStruct.h"
#   o "alnSetStruct.h"
#   o "alnSeqDefaults.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>  // by alnSetStructure.h
#   - <string.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o fun-01 WatermanSmithAln:
'    - Perform a Waterman Smith alignment on input
'      sequences
'  o fun-02 addBestBaseScore:
'    - Adds a score and index to the kept scores list
'  o fun-03 printMatrixCig:
'    - Prints out a cigar for an single path in a
'      direction matrix
'  o fun-04 printAltWaterAlns:
'    - Prints out the best aligment and the saved
'       alterantive alignments  (best alignment for each
'       base) to a file
'  o fun-05 updateDirScoreWaterSingle:
'    - Picks the best score and direction for the current
'      base pairs being compared in a Waterman-Smith
'      alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef WATERMAN_H
#define WATERMAN_H

#include <string.h>

#include "generalAlnFun.h"
#include "alnStruct.h"
#include "alnMatrixStruct.h"

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
      // seqST->offsetUI (index 0) and the point to end
      // the alignment seqST->endAlnUI (index 0).
      // CURRENTLY THIS ONLY WORKS FOR FULL ALIGNMENTS
    struct alnSet *settings,// Settings for the alignment
    // *startI and *endI paramaters should be index 1
    char *prefixCStr  // Prefix for matrix scan output file
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-01: addBestBaseScore
   '  - Adds a score and index to the kept scores list
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printMatrixCig
   '  - Prints out a cigar for an single path in a
   '    direction matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: updateDirScoreWaterSingle
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


#endif
