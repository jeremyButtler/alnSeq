/*#########################################################
# Name alignmentsFun
# Use:
#  o Holds functions for doing pairwise alignments (Needleman/waterman)
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
#   o <stdint.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  - fun-01 NeedleManWunschAln:
'     o Perform a Needleman-Wunsch alignment on input sequences
'  - fun-02 updateDirAndScore:
'     o Picks the best score and direction for the current base pairs
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef NEEDLEMAN_H
#define NEEDLEMAN_H

#include "generalAlnFun.h"
#include "alnStruct.h"
#include "alnMatrixStruct.h"

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: NeedlemanAln
   '  - Perform a Needleman-Wunsch alignment on input sequences
   '  o fun-01 sec-1: Variable declerations
   '  o fun-01 sec-2: Allocate memory for alignment
   '  o fun-01 sec-3: Fill in the initial negatives for the reference
   '  o fun-01 sec-4: Fill the matrix with scores
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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

#endif
