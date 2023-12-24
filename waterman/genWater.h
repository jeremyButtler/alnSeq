/*######################################################\
# Name: generalWater
#  - General functions for Waterman and memory efficent
#    Waterman alignments
# Libraries:
#  - "../general/genAln.h"
#  o "../general/genMath.h"
# C Standard Libraries:
\######################################################*/

#ifndef GENERAL_WATER_H
#define GENERAL_WATER_H

#include "../general/genAln.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o fun-01 waterMaxScore:
'    - Finds the max score and direction for a waterman
'      alignment
'  o fun-02 waterMaxEndRowScore:
'    - Finds the max score and direction for the last
'      score in a row for a waterman alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Fun-01: waterMaxScore
|   - Maximizes the score for a single waterman alignment
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
#define waterMaxScore(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,     /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
   insDir,      /*Insertion direction*/\
   nextSnpScore,/*Gets score to use for next snp*/\
   delScore,    /*Score for an deletion*/\
   alnSetPtr    /*Pointer to alnSet with settings*/\
){ /*waterMaxScore*/\
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
    charMaxScore(\
       (scoreOn),\
       (dirOn),\
       (macroSnpScoreL),\
       (macroInsScoreL),\
       (delScore),\
       (alnSetPtr)->bestDirC\
    ); /*Find the best direction*/\
    \
    /*Find if keeping score of if alignment stops*/\
    keepDirUL = -((scoreOn) > 0);\
    (dirOn) &= keepDirUL;\
    (scoreOn) &= keepDirUL;\
    \
    (delScore) =\
         (scoreOn)\
       + ((gapDiff) & (-(dirOn) != defMvSnp))\
       + (alnSetPtr)->gapOpenC;\
} /*waterMaxScore*/

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
#define waterMaxEndRowScore(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,    /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
   insDir,      /*Insertion direction*/\
   nextSnpScore,/*Gets score to use for next snp*/\
   delScore,    /*Score for an deletion*/\
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
    charMaxScore(\
       (scoreOn),\
       (dirOn),\
       (macroSnpScoreL),\
       (macroInsScoreL),\
       (delScore),\
       (alnSetPtr)->bestDirC\
    ); /*Find the best direction*/\
    \
    /*Find if keeping score of if alignment stops*/\
    keepDirUL = -((scoreOn) > 0);\
    (dirOn) &= keepDirUL;\
    (scoreOn) &= keepDirUL;\
} /*waterMaxEndRowScore*/

#endif