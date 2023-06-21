/*######################################################################
# Name alignmentsFun
# Use:
#  o Holds functions for doing pairwise alignments (Needleman/waterman)
# Includes:
#   - "alnSeqDefaults.h"
#   - "cStrToNumberFun.h"
#   - "twoBitArrays.h"
# C Standard libraries:
#   - <stdio.h>
#   o <stdlib.h>
#   o <stdint.h>
######################################################################*/

#ifndef ALIGNMENTSFUN_H
#define ALIGNMENTSFUN_H

#include "alnSeqDefaults.h"
#include "cStrToNumberFun.h"
#include "twoBitArrays.h"

#include <stdio.h>

#define defClearNonAlph (1 | 2 | 4 | 8 | 16) // clear 64 bit and case
#define defToUper (1 | 2 | 4 | 8 | 16 | 64)
    // Clear 32nd bit (marks lower case)

#define defMoveStop 0    // Do not move
#define defMoveLeft 1    // Move left (deletion) in alignment matrix
#define defMoveUp 2      // Move up (insertion) in alignment matrix
#define defMoveDiagnol 3 // Move on a diagnol (snp/match) in alignment

#define defDelFlag 1    // deletion
#define defInsFlag 2    // insertion
#define defBaseFlag 4   // match or snp
#define defSoftQueryFlag 8 // Softmask a query base
#define defSoftRefFlag 16  // Softmask a reference base

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o st-01 alnSet:
'     o Holds settings for my alignment program
'  o fun-01 initAlnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
'  o fun-02 freeAlnSet:
'    o Frees and alnSet (alignment settings) structure
'  o fun-03 setBasePairScore:
'    - Changes SNP/Match penalty for one query/reference
'      combination
'  - fun-04 readInScoreFile
'     o Reads in a file of scores for a scoring matrix
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| ST-01: alnSet
| Use: Holds settings for my alignment program
\---------------------------------------------------------------------*/
typedef struct alnSet
{ /*alnSet*/
   // Line wrap for printing out an alignment
   unsigned short lineWrapUS;
   unsigned short lenFileNameUS;

   //kmer mapping variables
   char diagnolPriorityC; // 0 favor snps; 1 kinda; 2 do not
   char topPriorityC;     // 0 favor insertions; 1 kinda; 2 do not
   char leftPriorityC;    // 0 favor deletions; 1 kinda; 2 do not

   // Needleman-Wunsch / Waterman Smith variables
   char useNeedleBl;
   char multiBaseWaterBl; // Keep more than best aligment
   char refQueryScanBl;   // Keep best score for ref/query
   char matrixScanBl;  // Fur doing a matrix scan
     // If set to 1: Recored a best score for each base
     //   in the reference and query in a Smith Waterman
     // alignment
   int16_t snpPenaltyC[26][26];   // Penalty for mismatches (in matrix)
     // Size is due to wanting a look up table that can handle
     // anonymous bases. Most cells will be set to 0.
     // value = snpPenaltyC[(uint8_t) (base1 & defClearNonAlph) - 1 ]
     //                    [(uint8_t) (base2 & defClearNonAlph) - 1 ]
   int32_t gapStartPenaltyI;     // Penalty for starting an indel
   int32_t gapExtendPenaltyI;    // Penalty for extending an indel
   uint32_t minScoreUI;        // Minimum score needed to keep alignment
   uint32_t minBasesUI;
}alnSet;

/*---------------------------------------------------------------------\
| Output: Modifies: alnSetST to have default alignment settings values
\---------------------------------------------------------------------*/
void initAlnSet(
    struct alnSet *alnSetST // Alinment settings structure to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: initAlnSet
   '  - Set values in altSet (alingment settings) structure to defaults
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: one score in an snp/match scoring matrix
\---------------------------------------------------------------------*/
void setBasePairScore(
    const char *queryBaseC, // Query base to change score for
    const char *refBaseC,   // Reference base to change score for
    int16_t newScoreC,      // New value for [query][ref] combination
    struct alnSet *alnSetST // structure with scoring matrix to change
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: setBasePairScore
   '  - Changes SNP/Match penalty for one query/reference combination
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Frees the alnSet structure (does not set pointer to 0)
\---------------------------------------------------------------------*/
void freeAlnSet(
    struct alnSet *alnSetST,  // Alignment settings structure to free
    char stackBl              // 1: alnSetSt on stack; 0: on heap
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

unsigned long readInScoreFile(
    struct alnSet *alnSetST,  // structure with scoring matrix to change
    FILE *scoreFILE           // File of scores for a scoring matrix
);  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-04 TOC: readInScoreFile
    '  - Reads in a file of scores for a scoring matrix
    '  o fun-04 sec-1: Variable declerations and buffer set up
    '  o fun-04 sec-2: Read in line and check if comment
    '  o fun-04 sec-3: Get score, update matrix, & move to next line
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
