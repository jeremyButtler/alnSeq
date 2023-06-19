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
'  - st-01 alnSet:
'     o Holds settings for my alignment program
'  - fun-03 cnvtAlnErrToSeq:
'     o Uses an array of error types and a c-string with the a sequence
'       to make one part of an alignment
'  - fun-04 WatermanSmithAln:
'     o Perform a Waterman Smith alignment on input sequences
'  - fun-05 NeedleManWunschAln:
'     o Perform a Needleman-Wunsch alignment on input sequences
'  - fun-07 checkIfBasesMatch
'     o Check if two bases are the same (includes anonymous bases)
'  - fun-12 makeAlnSet:
'     o Makes & initalizes an alnSet structer on the heap
'  - fun-13 getBasePairScore:
'     o Get the score for a pair of bases from an alignment structure
'  - fun-14 initAlnSet:
'     o Set values in altSet (alingment settings) structure to defaults
'  - fun-15 setAlnSetSnpPenalty:
'     o Changes SNP/Match penalty for one query/reference combination
'  - fun-16 freeAlnSet:
'     o Frees and alnSet (alignment settings) structure
'  - fun-17 cnvtAlnErrAryToLetter:
'     o Converts an alignment error array from my Needleman-Wunsch
'       alignment into an array of letters (I = insertion, D = deletion,
'       = = match, X = snp) [These codes are from the eqx cigar entry]
'  - fun-18 readInScoreFile
'     o Reads in a file of scores for a scoring matrix
'  - fun-19 getIndelScore:
'     o Gets an indel score for the current cell
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
| Output: Heap alloacted C-string with alignment for the input sequence
\---------------------------------------------------------------------*/
char * cnvtAlnErrToSeq(
    char *seqCStr,        // c-string with query sequence
    int32_t seqStartI,    // Were alignment starts on query (index 1)
    char queryBl,         // 1: input sequence is query, 0: is reference
    uint8_t *alnErrUCAry, // Holds error types (ends with 0)
    uint32_t lenErrAryUI  // length of alnErrUCAry (for sequence array)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: cnvtAlnQueryErrToSeq
   '  - Uses an array of error types and a c-string with the a sequence
   '    to make one part of an alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Returns: Initalized alnSet structer or 0 for memory error
| Note: I prefer using stack, so this is more here for someone else
\---------------------------------------------------------------------*/
struct alnSet * makeAlnSetST(
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: makeAlnSet
   '  - Makes & initalizes an alnSet structer on the heap
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: alnSetST to have default alignment settings values
\---------------------------------------------------------------------*/
void initAlnSet(
    struct alnSet *alnSetST // Alinment settings structure to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: initAlnSet
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
   ' Fun-15 TOC: Sec-1 Sub-1: setBasePairScore
   '  - Changes SNP/Match penalty for one query/reference combination
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Frees the alnSet structure (does not set pointer to 0)
\---------------------------------------------------------------------*/
void freeAlnSet(
    struct alnSet *alnSetST,  // Alignment settings structure to free
    char stackBl              // 1: alnSetSt on stack; 0: on heap
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: alnErrUCAry to have letters instead of codes
\---------------------------------------------------------------------*/
void cnvtAlnErrAryToLetter(
    char *refSeqCStr,      // Reference sequence for detecting matches
                           // Input 0 to ignore matches
    char *querySeqCStr,    // query sequence for detecting matches
                           // Input 0 to ignore matches
    uint8_t *alnErrUCAry   // Array array from Needleman (fun-05) to 
                          // convert to legible characters
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-17 TOC: Sec-1 Sub-1: cnvtAlnErrAryToLetter
   '  - Converts an alignment error array from my Needleman-Wunsch
   '    alignment into an array of letters (I = insertion, D = deletion,
   '    = = match, X = snp) [These codes are from the eqx cigar entry]
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

unsigned long readInScoreFile(
    struct alnSet *alnSetST,  // structure with scoring matrix to change
    FILE *scoreFILE           // File of scores for a scoring matrix
);  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-18 TOC: readInScoreFile
    '  - Reads in a file of scores for a scoring matrix
    '  o fun-18 sec-1: Variable declerations and buffer set up
    '  o fun-18 sec-2: Read in line and check if comment
    '  o fun-18 sec-3: Get score, update matrix, & move to next line
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
