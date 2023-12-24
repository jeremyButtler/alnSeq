/*#######################################################\
# Name: generalQryRefScan
#   - Holds functions to print and manage query refernce
#     scans 
# Libraries:
#   - "../general/genScan.h"           (No .c file)
#   o "../general/genMath.h"           (No .c file)
#   o "../general/alnMatrixStruct.h"   (No .c file)
#   o "../general/twoBitArrays.h"      (No .c file)
#   o "../general/dataTypeShortHand.h" (No .c file)
# C Standard Libraries:
#   o <stdio.h>
\#######################################################*/


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  - Functions and variables for working with the query
'    reference scan Waterman alignments.
'  o header:
'    - definitions and includes
'  o fun-01 waterMatrixScanMaxScore:
'    - Maximizes the score for a single base pair in an
'      waterman query reference scan
'  o fun-02 waterMatrixScanMaxEndRow:
'    - Maximizes the score for a single waterman alignment
'    - This for the end of a row, were you do not want the
'      next snp or deletion score
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*------------------------------------------------------\
| Header: 
|   - Definitions and includes
\------------------------------------------------------*/

#ifndef GEN_WATER_SCAN_H
#define GEN_WATER_SCAN_H

#include "../general/genScan.h"

/*-------------------------------------------------------\
| Fun-01: waterScanMaxScore
|   - Maximizes the score for a single base pair in an
|     waterman query reference scan
| Input;
|   - refBase;
|     o Reference base to find score for
|   - qryBase;
|     o Query base to find score for
|   - gapDiff:
|     o Gap Exentsion score - gap open score.
|     o This is used to find the indel scores
|   - scoreOn:
|     o Score to update/score for an insertion
|   - dirOn;
|     o Direction to update/direction did for last ins
|   - nextSnpScore:
|     o Score to use in finding the snp score
|   - delScore:
|     o Score for making a deletion (is updated)
|   - alnSetPtr:
|     o Pointer to alnSet structure with the settings
|       for the alignment
| Output:
|   - Modifies:
|     o socreOn to hold the maximum score
|     o dirOn to hold the direction of the max score
|     o nextSnpScore to hold the score used to find the
|       next snp score
|     o delScore to hold the score for the next deletion
\-------------------------------------------------------*/
#define waterMatrixScanMaxScore(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,    /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
   insDir,      /*Previous ins direction*/\
   index,       /*Start of alignment for kept dir*/\
   nextSnpScore,/*Gets score to use for next snp*/\
   delScore,    /*Score for an deletion*/\
   snpPos,      /*Index for snp*/\
   insPos,      /*Indes for ins*/\
   delPos,      /*Index for del*/\
   curIndex,    /*Current index (for stops)*/\
   alnSetPtr    /*Pointer to alnSet with settings*/\
){ /*waterScanMaxScore*/\
   long macroSnpScoreL =\
        (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
    \
    long macroInsScoreL =\
         (scoreOn)\
       + ((gapDiff) & (-((insDir) != defMvSnp)))\
       + (alnSetPtr)->gapOpenC;\
    \
    ulong keepDirUL = 0;\
    \
    (nextSnpScore) = (scoreOn);/*Score to find next snp*/\
    \
    scanMaxScore(\
       (scoreOn),\
       (dirOn),\
       (index),\
       (macroSnpScoreL),\
       (macroInsScoreL),\
       (delScore),\
       (snpPos),\
       (insPos),\
       (delPos),\
       (alnSetPtr)->bestDirC\
    ); /*Find the best direction*/\
    \
    /*Find if keeping score of if alignment stops*/\
    keepDirUL = -((scoreOn) > 0);\
    (dirOn) &= keepDirUL;\
    (scoreOn) &= keepDirUL;\
    (index) =\
         ((index) & keepDirUL)\
       + ((curIndex) & (~keepDirUL));\
    \
    (delScore) =\
         (scoreOn)\
       + ((gapDiff) & (-(dirOn) != defMvSnp))\
       + (alnSetPtr)->gapOpenC;\
} /*waterScanMaxScore*/

/*-------------------------------------------------------\
| Fun-02: waterMaxEndRowScore
|   - Maximizes the score for a single waterman alignment
|   - This for the end of a row, were you do not want the
|     next snp or deletion score
| Input;
|   - refBase;
|     o Reference base to find score for
|   - qryBase;
|     o Query base to find score for
|   - gapDiff:
|     o Gap Exentsion score - gap open score.
|     o This is used to find the indel scores
|   - scoreOn:
|     o Score to update/score for an insertion
|   - dirOn;
|     o Direction to update/direction did for last ins
|   - nextSnpScore:
|     o Score to use in finding the snp score
|   - alnSetPtr:
|     o Pointer to alnSet structure with the settings
|       for the alignment
| Output:
|   - Modifies:
|     o socreOn to hold the maximum score
|     o dirOn to hold the direction of the max score
\-------------------------------------------------------*/
#define waterMatrixScanMaxEndRow(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,    /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
   insDir,      /*Direction of last insertion*/\
   index,       /*Start of alignment for kept dir*/\
   nextSnpScore,/*Gets score to use for next snp*/\
   delScore,    /*Score for an deletion*/\
   snpPos,      /*Index for snp*/\
   insPos,      /*Indes for ins*/\
   delPos,      /*Index for del*/\
   curIndex,    /*Current index (for stops)*/\
   alnSetPtr    /*Pointer to alnSet with settings*/\
){ /*waterMaxEndRowScore*/\
   long macroSnpScoreL =\
      (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
    \
    long macroInsScoreL =\
         (scoreOn)\
       + ((gapDiff) & (-((insDir) != defMvSnp)))\
       + (alnSetPtr)->gapOpenC;\
    \
    ulong keepDirUL = 0;\
    \
    scanMaxScore(\
       (scoreOn),\
       (dirOn),\
       (index),\
       (macroSnpScoreL),\
       (macroInsScoreL),\
       (delScore),\
       (snpPos),\
       (insPos),\
       (delPos),\
       (alnSetPtr)->bestDirC\
    ); /*Find the best direction*/\
    \
    /*Find if keeping score of if alignment stops*/\
    keepDirUL = -((scoreOn) > 0);\
    (dirOn) &= keepDirUL;\
    (scoreOn) &= keepDirUL;\
    (index) =\
         ((index) & keepDirUL)\
       + ((curIndex) & (~keepDirUL));\
} /*waterMaxEndRowScore*/

#endif
