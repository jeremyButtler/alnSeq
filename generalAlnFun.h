/*########################################################
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
#   o <stdio.h>  // Used by alnSetStructures.h
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o fun-01 getIndelScore:
'    - Gets an indel score for the current cell
'  o fun-02 getBasePairScore:
'    - Get the score for a pair of bases from an alignment
'  o fun-03 checkIfBasesMatch:
'    - Are two bases are same? (includes anonymous bases)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef GENERALLALNFUN_H
#define GENERALLALNFUN_H

#include "alnSetStruct.h" // Default settings

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o the score for an indel
\--------------------------------------------------------*/
long getIndelScore(
    struct twoBitAry *lastDir, // Has last direction
    struct alnSet *alnSetST,
      // Has gap open & extension penalties
    long *lastBaseL      // Has score of the last base
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: getIndelScore
   '  - Gets an indel score for the current cell
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o score of a single pair of bases
\--------------------------------------------------------*/
int16_t getBasePairScore(
    const char *queryBaseC, // Query base to score
    const char *refBaseC,   // Reference base to score
    struct alnSet *alnSetST // has scoring matrix
);/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ' Fun-02 TOC: Sec-1 Sub-1: getBasePairScore
  '  - Get the score for a pair of bases from an alignment
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o 1: if bases were a match
|    o 0 if bases do not mach
\--------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC,// Query base to compare to reference
    char *refBaseC   // Reference base to compare to query
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Are two bases are same? (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
