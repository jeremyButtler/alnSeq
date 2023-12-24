/*######################################################\
# Name: generalWater
#  - General functions for Waterman and memory efficent
#    Waterman alignments
# Libraries:
#  - "../general/genAln.h"
#  o TODO: FILL OUT
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdio.h>
#  o TODO: FILL OUT
\######################################################*/

#ifndef GENERAL_NEEDLE_H
#define GENERAL_NEEDLE_H

#include "../general/genAln.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o fun-01 needleMaxScore:
'    - Finds the max score and direction for an Needleman
'      alignment
'  o fun-02 needleMaxEndRowScore:
'    - Finds the max score and direction for the last
'      score in a row for an Needleman alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Fun-01: needleMaxScore
|   - Maximizes the score for a single Needleman alignment
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
#define needleMaxScore(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,     /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
   insDir,      /*Insertion direction*/\
   nextSnpScore,/*Gets score to use for next snp*/\
   delScore,    /*Score for an deletion*/\
   alnSetPtr    /*Pointer to alnSet with settings*/\
){ /*needleMaxScore*/\
   long macroSnpScoreL =\
        (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
    \
    long macroInsScoreL =\
         (scoreOn)\
       + ((gapDiff) & (-((insDir) != defMvSnp)))\
       + (alnSetPtr)->gapOpenC;\
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
    (delScore) =\
         (long) (scoreOn)\
       + ((gapDiff) & (-(dirOn) != defMvSnp))\
       + (alnSetPtr)->gapOpenC;\
} /*needleMaxScore*/

/*-------------------------------------------------------\
| Fun-02: needleMaxEndRowScore
|   - Maximizes the score for a single Needleman alignment
|   - This for the end of a row, were you do not want the
|     next snp or deletion score
| Input;
|   - refBase;
|     o Reference base to find score for
|   - qryBase;
|     o Query base to find score for
|   - gapDiff:j
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
#define needleMaxEndRowScore(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,    /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
   insDir,      /*Insertion direction*/\
   nextSnpScore,/*Gets score to use for next snp*/\
   delScore,    /*Score for an deletion*/\
   alnSetPtr    /*Pointer to alnSet with settings*/\
){ /*needleMaxEndRowScore*/\
   long macroSnpScoreL =\
      (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
    \
    long macroInsScoreL =\
         (scoreOn)\
       + ((gapDiff) & (-((insDir) != defMvSnp)))\
       + (alnSetPtr)->gapOpenC;\
    \
    charMaxScore(\
       (scoreOn),\
       (dirOn),\
       (macroSnpScoreL),\
       (macroInsScoreL),\
       (delScore),\
       (alnSetPtr)->bestDirC\
    ); /*Find the best direction*/\
} /*needleMaxEndRowScore*/

#endif