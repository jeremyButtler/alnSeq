/*#######################################################\
# Name: generalQryRefScanNoGap
#   - Holds functions to print and manage query refernce
#     scans without using gap extension penalties
# Libraries:
#   - "../general/genWaterScanNoGap.h" (No .c file)
#   o "../general/alnMatrixStruct.h"   (No .c file)
#   o "../general/genMath.h"           (No .c file)
#   o "../general/twoBitArrays.h"      (No .c file)
#   o "../general/dataTypeShortHand.h" (No .c file)
# C Standard Libraries:
#   o <stdio.h>
\#######################################################*/


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  - Functions and variables for working with the query
'    reference scan Waterman alignments. This does not
'    use gap extension penalties
'  o header:
'    - definitions and includes
'  o fun-01 waterMatrixScanMaxScoreNoGap:
'    - Maximizes the score for a single base pair in an
'      waterman query reference scan. This does not use
'      gap extension penalties
'  o fun-02 waterMatrixScanMaxEndRowNoGap:
'    - Maximizes the score for a single waterman alignment
'    - This for the end of a row, were you do not want the
'      next snp or deletion score. This does not use gap
'      extension penalties
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*------------------------------------------------------\
| Header: 
|   - Definitions and includes
\------------------------------------------------------*/

#ifndef GEN_WATER_SCAN_NO_GAP_H
#define GEN_WATER_SCAN_NO_GAP_H

#include "../general/genScan.h"

/*-------------------------------------------------------\
| Fun-01: waterScanMaxScoreNoGap
|   - Maximizes the score for a single base pair in an
|     waterman query reference scan. This does not use
|     gap extension penalties
| Input;
|   - refBase;
|     o Reference base to find score for
|   - qryBase;
|     o Query base to find score for
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
#define waterMatrixScanMaxScoreNoGap(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
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
){ /*waterScanMaxScoreNoGap*/\
   long macroSnpScoreL =\
        (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
    \
    long macroInsScoreL =\
         (scoreOn) + (alnSetPtr)->gapOpenC;\
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
         (scoreOn) + (alnSetPtr)->gapOpenC;\
} /*waterScanMaxScoreNoGap*/

/*-------------------------------------------------------\
| Fun-02: waterMaxEndRowScoreNoGap
|   - Maximizes the score for a single waterman alignment
|     this does not use gap extension penalties
|   - This for the end of a row, were you do not want the
|     next snp or deletion score
| Input;
|   - refBase;
|     o Reference base to find score for
|   - qryBase;
|     o Query base to find score for
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
#define waterMatrixScanMaxEndRowNoGap(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
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
){ /*waterMaxEndRowScoreNoGap*/\
   long macroSnpScoreL =\
      (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
    \
    long macroInsScoreL =\
         (scoreOn) + (alnSetPtr)->gapOpenC;\
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
} /*waterMaxEndRowScoreNoGap*/

#endif
