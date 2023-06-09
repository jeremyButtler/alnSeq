/*#########################################################
# Name generalAlnFun
# Use:
#  o Holds general functions used in my Needleman Wunsch
#    or Waterman Smith alignment.
# Libraries:
#   - "alnSetStruct.h"
#   o "alignmentSettings.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>  // Used by alnSetStructure.h
#########################################################*/

#include "generalAlnFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 getIndelScore:
'    - Gets an indel score for the current cell
'  o fun-02 getBasePairScore:
'    - Get the score for a pair of bases from an alignment
'  o fun-03 cnvtAlnErrToSeq:
'    - Uses an array of error types and a c-string with
'      the a sequence to make one part of an alignment
'  o fun-04 TOC: Sec-1 Sub-1: cnvtAlnErrAryToLetter
'    - Converts an alignment error array from my alignment
'      algorithms into an array of letters
'      (I = insertion, D = deletion, = = match, X = snp)
'  o fun-05 checkIfBasesMatch:
'    - Are two bases are same? (includes anonymous bases)
'  o fun-06 getAlnAry:
'    - Builds an anlignment array for the input direction
'      matrix
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o the score for an indel
\--------------------------------------------------------*/
long getIndelScore(
    uint8_t *lastDirUC,  // Cell gettign last score from
    uint8_t *lastBitUC,  // two bit element on in lastDirUC
    struct alnSet *alnSetST,
      // Has gap open & extension penalties
    long *lastBaseL      // Has score of the last base
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: getIndelScore
   '  - Gets an indel score for the current cell
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Check if this is the first indel or not
   switch(getTwoBitAryElm(lastDirUC, lastBitUC))
   { // Switch; check if this is the first indel
     case defMoveStop:     // At the end of the matrix
     case defMoveDiagnol:  // Top base was a SNP
         return *lastBaseL + alnSetST->gapStartPenaltyI;

     case defMoveUp:   // Top base was an insertion
     case defMoveLeft: // Top base was an deletion
         return *lastBaseL + alnSetST->gapExtendPenaltyI;
   } // Switch; check if this is the first indel

   return 0; // Somthing went wrong
} // getIndelScore

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o score of a single pair of bases
\--------------------------------------------------------*/
int16_t getBasePairScore(
    const char *queryBaseC, // Query base to score
    const char *refBaseC,   // Reference base to score
    struct alnSet *alnSetST // has scoring matrix
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ' Fun-02 TOC: Sec-1 Sub-1: getBasePairScore
  '  - Get the score for a pair of bases from an alignment
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   return
       alnSetST->snpPenaltyC
           [(uint8_t) (*queryBaseC & defClearNonAlph) - 1]
           [(uint8_t) (*refBaseC & defClearNonAlph) - 1];
} // getBasePairScore

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o 1: if bases were a match
|    o 0 if bases do not mach
\--------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC,// Query base to compare to reference
    char *refBaseC   // Reference base to compare to query
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Are two bases are same? (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
   // The switch should default to a look up table &
   // will be more clear
   switch(*queryBaseC & defToUper)
   { // Switch: Check if bases are the same
       case 'A':
       // Case: Query is an A
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 1;
               case 'W': return 1;
               case 'M': return 1;
               case 'R': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an A

       case 'T':
       // Case: Query is an T
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'T': return 1;
               case 'U': return 1;
               case 'W': return 1;
               case 'K': return 1;
               case 'B': return 1;
               case 'Y': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an T

       case 'U':
       // Case: Query is an U
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'T': return 1;
               case 'U': return 1;
               case 'W': return 1;
               case 'K': return 1;
               case 'B': return 1;
               case 'Y': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an U

       case 'C':
       // Case: Query is an C
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 1;
               case 'S': return 1;
               case 'M': return 1;
               case 'Y': return 1;
               case 'B': return 1;
               case 'H': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an C

       case 'G':
       // Case: Query is an G
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'G': return 1;
               case 'S': return 1;
               case 'K': return 1;
               case 'R': return 1;
               case 'B': return 1;
               case 'D': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an G

       case 'W':
       // Case: Query is an W
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 0;
               case 'G': return 0;
               case 'S': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an W

       case 'S':
       // Case: Query is an S
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'W': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an S

       case 'M':
       // Case: Query is an M
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'G': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'K': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an M

       case 'K':
       // Case: Query is an K
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               case 'C': return 0;
               case 'M': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an K

       case 'R':
       // Case: Query is an R
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'Y': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an R

       case 'Y':
       // Case: Query is an Y
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               case 'G': return 0;
               case 'R': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an Y

       case 'B':
       // Case: Query is an B
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an B

       case 'D':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an C

       case 'H':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'G': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an H

       case 'V':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'T': return 0;
               case 'U': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an V

       case 'N': return 1;
       case 'X': return 1;
   } // Switch: Check if bases are the same

   return 0; // not a base
} // checkIfBasesMatch
