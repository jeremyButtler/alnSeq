/*########################################################
# Name: generalScoreHirsch
# Use:
#  - Holds functions for doing scoring in a scalar
#    Hirschberg alignment
# Libraries:
#  - "genHirsch.h"                    (No .c File)
#  o "../general/alnStruct.h"         (No .c File)
#  o "../general/seqStruct.h"         (No .c File)
#  o "../general/alnSetStruct.h"      (No .c File)
#  o "../general/alnSeqDefaults.h"    (No .c File)
#  o "../general/strToNum.h"          (No .c File)
#  o "../general/dataTypeShortHand.h" (No .c File)
#  o "../general/alnMatrixStruct.h"   (No .c File)
#  o "../general/genAln.h"            (No .c File)
#  o "../general/genMath.h"           (No .c File)
# C Standard Libraries:
#  o <time.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <string.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - Includes and definitions
'  o fun-01 hirschScore:
'    - Gets a single score for a single base pair in an
'      Hirschberg alignment
'  o fun-02 hirschScoreRowEnd:
'    - Finds the last score for the last base pair in a
'      row for a Hirschberg alignment
'  o fun-03 scoreHirschFor:
'    - Does a single round of scoring for a hirschberg
'      alignment (forward direction)
'  o fun-04 scoreHirschRev:
'    - Does a single round of scoring for a hirschberg
'      alignment (reverse direction)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Includes and definitions
\-------------------------------------------------------*/

#ifndef GENERAL_SCORE_HIRSCH_H 
#define GENERAL_SCORE_HIRSCH_H 

#include "genHirsch.h"

/*-------------------------------------------------------\
| Fun-01: hirschScore
|   - Gets a single score for a single base pair in an
|     Hirschberg alignment
| Input:
|   - refBase:
|     o Reference base to get the score for (character)
|   - qryBase:
|     o Query base to get the score for (character)
|   - gapDiff:
|     o The gap extension penatly - gap starting penatly.
|       This is used to find the insertion and deltion
|       scores
|   - scoreOn:
|     o Score to updated and used to find insertions
|   - dirOn:
|     o Boolean with the direction on.
|       -1: for a gap
|        0: for not a gap
|     o Is updated to hold the new direction
|   - nextSnpScore:
|     o Has the score for finding an snp.
|     o This is updated to hold the next snps score
|   - delScore:
|     o Score for a deletion
|     o This is updated to hold the next deletion
|   - alnSetPtr:
|     o Pointer to an alnSet structer with alnSetPtr to
|       use in the alignment
| Output:
|   - Modifies:
|     o scoreOn to hold the next maximum score
|     o dirOn to hold the next direction
|     o nextSnpScore to hold the next score needed to
|       calculate the next snp
|     o delScore to hold the next deletion score
\-------------------------------------------------------*/
#define hirschScore(\
   refBase,  /*Reference base*/\
   qryBase,  /*Query sequence*/\
   gapDiff,  /*gapExtend - gapOpen; used for gap score*/\
   scoreOn,  /*Score on and holds the new maximum*/\
   dirOn,    /*Current direction & holds next direction*/\
   nextSnpScore, /*Holds score used to get snp score*/\
   delScore,  /*Holds the deletion score*/\
   alnSetPtr  /*Holds the alnSetPtr for the alignment*/\
){ /*hirschScore*/\
   long macroSnpScoreL =\
        (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
   \
   /*Get the insertion score*/\
   long macroInsScoreL =\
        (scoreOn)\
      + ((gapDiff) & (dirOn)) + (alnSetPtr)->gapOpenC;\
      /* Logic (ins/del): (gapDiff & dir) + gapOpen
      `   o gapDiff = gapExtend - gapOpen
      `   o dirOn is -1 or 0 (Nothing else is valid)
      `       - if -1; gapDiff & -1 = gapDiff
      `       - if 0; gapDiff & 0 = 0
      `   o gapDiff + gapOpen = gapExtend
      `   o 0 + gapOpen = gapOpen
      */\
   \
   /*Get the score needed to find the next snp/match. This
   ` is overwiten when I find the max score.
   */\
   (nextSnpScore) = (scoreOn);\
   \
   maxGapScore(\
      (scoreOn),\
      (dirOn),\
      macroSnpScoreL,           /*snp/match score*/\
      macroInsScoreL,           /*insertion score*/\
      (delScore),               /*deletion score*/\
      (alnSetPtr)->bestDirC\
   );\
   \
   (delScore) =\
        (scoreOn)\
      + ((gapDiff) & (dirOn)) + (alnSetPtr)->gapOpenC;\
} /*hirschScore*/

/*-------------------------------------------------------\
| Fun-02: hirschScoreEnd
|   - Finds the last score for the last base pair in a
|     row for a Hirschberg alignment
| Input:
|   - refBase:
|     o Reference base to get the score for (character)
|   - qryBase:
|     o Query base to get the score for (character)
|   - gapDiff:
|     o The gap extension penatly - gap starting penatly.
|       This is used to find the insertion and deltion
|       scores
|   - scoreOn:
|     o Score to update and is used to find insertions
|   - dirOn:
|     o Boolean with the direction on.
|       -1: for a gap
|        0: for not a gap
|     o Is updated to hold the new direction
|   - nextSnpScore:
|     o Has the score for finding an snp.
|   - delScore:
|     o Score for a deletion
|   - alnSetPtr:
|     o Pointer to an alnSet structer with alnSetPtr to
|       use in the alignment
| Output:
|   - Modifies:
|     o scoreOn to hold the next maximum score
|     o dirOn to hold the next direction
\-------------------------------------------------------*/
#define hirschScoreRowEnd(\
   refBase,  /*Reference base*/\
   qryBase,  /*Query sequence*/\
   gapDiff,  /*gapExtend - gapOpen; used for gap score*/\
   scoreOn,  /*Holds new maxium score*/\
   dirOn,    /*Current direction (updated)*/\
   nextSnpScore, /*Holds score used to get snp score*/\
   delScore,  /*Holds the deletion score*/\
   alnSetPtr  /*Holds the alnSetPtr for the alignment*/\
){ /*hirschScoreRowEnd*/\
   long macroSnpScoreL =\
        (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
   \
   /*Get the insertion score*/\
   long macroInsScoreL =\
        (scoreOn)\
      + ((gapDiff) & (dirOn)) + (alnSetPtr)->gapOpenC;\
      /* Logic (ins/del): (gapDiff & dir) + gapOpen
      `   o gapDiff = gapExtend - gapOpen
      `   o dirOn is -1 or 0 (Nothing else is valid)
      `       - if -1; gapDiff & -1 = gapDiff
      `       - if 0; gapDiff & 0 = 0
      `   o gapDiff + gapOpen = gapExtend
      `   o 0 + gapOpen = gapOpen
      */\
   maxGapScore(\
      (scoreOn),\
      (dirOn),\
      macroSnpScoreL,           /*snp/match score*/\
      macroInsScoreL,           /*insertion score*/\
      (delScore),               /*deletion score*/\
      (alnSetPtr)->bestDirC\
   );\
} /*hirschScoreRowEnd*/

/*-------------------------------------------------------\
| Fun-03: scoreHirschFor
|  - Scores the input query and reference sequence
| Input:
|  - refSeqStr:
|    o C-strint with the reference sequence to align
|  - refStart:
|    o Position of the first base to align in the
|      reference sequence (is index 0)
|  - refLen:
|    o Number of bases to align in the reference sequence
|      (is index 1)
|  - qrySeqStr:
|    o C-string with the query sequence to align
|  - qryStart:
|    o Position of the first base to align in the query
|      sequence (is index 0)
|  - qryLen:
|    o Number of bases to align in the query sequence
|      (is index 1)
|  - scoreAry:
|    o Array to fill with scores
|  - dirAry:
|    o array to use in finding directions
|  - alnSetPtr:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o scoreAry to hold the last row of scores in a
|      Needleman Wunsch
|    o dirAry to hold the last row of directions in a
|      Needleman Wunsch
|  - Returns:
|    o The gap column score
\-------------------------------------------------------*/
#define scoreHirschFor(\
  refSeqStr,   /*Reference sequence*/\
  refStart,    /*index 0 starting ref base*/\
  refLen,      /*index 1 Length of ref*/\
  \
  qrySeqStr,   /*Query sequence*/\
  qryStart,    /*Index 0 Starting query base*/\
  qryLen,      /*index 1 length of query*/\
  \
  scoreAry,    /*Array of scores to fill*/\
  dirAry,      /*direction row, for gap extend*/\
  alnSetPtr    /*setttings to use*/\
)({/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: scoreHirschFor
   '  - Does the scoring step for a hirschberg alignment
   '    (forward direction)
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-03 sec-03:
   '    - Score till on the last row
   '  o fun-03 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   short gapDiffS =\
      (alnSetPtr)->gapExtendC - (alnSetPtr)->gapOpenC;\
      /*Used to find the gap score (branchless method)*/\
   \
   ulong refBaseMac = 0;\
   ulong qryBaseMac = 0;\
   long macroGapColL = 0;\
   long macro_score_for_next_snp_L = 0;\
   long macroDelScoreL = 0;\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Set up the first row (indel row) of scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   macroGapColL = 0;\
   (scoreAry)[refStart] = (alnSetPtr)->gapOpenC;\
   (dirAry)[refStart] = -1;\
   \
   for(\
      refBaseMac = (refStart) + 1;\
      refBaseMac < (refStart) + (refLen);\
      ++refBaseMac\
   ){ /*Loop:Set the initial blank scores*/\
      (dirAry)[refBaseMac] = -1;\
      \
      (scoreAry)[refBaseMac] =\
           (scoreAry)[refBaseMac - 1]\
         + (alnSetPtr)->gapExtendC;\
   } /*Loop:Set the initial blank scores*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Find the scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*Get the value I need to find the next snp's score*/\
   macro_score_for_next_snp_L = 0;\
   macroGapColL = (alnSetPtr)->gapOpenC;\
   macroDelScoreL = macroGapColL + alnSetPtr->gapExtendC;\
   \
   for(\
      qryBaseMac = (qryStart);\
      qryBaseMac < (qryStart) + (qryLen);\
      ++qryBaseMac\
   ){ /*Loop: score all query bases (rows)*/\
      \
      for(\
         refBaseMac = (refStart);\
         refBaseMac < (refStart) + (refLen) - 1;\
         ++refBaseMac\
      ){ /*Loop: Find the max scores for a single row*/\
        hirschScore(\
           (refSeqStr)[refBaseMac],\
           (qrySeqStr)[qryBaseMac],\
           gapDiffS,\
           (scoreAry)[refBaseMac],\
           (dirAry)[refBaseMac],\
           macro_score_for_next_snp_L,\
           macroDelScoreL,\
           (alnSetPtr)\
        ); /*Get max score for the current base pair*/\
      } /*Loop: Find the max scores for a single row*/\
      \
      hirschScoreRowEnd(\
         (refSeqStr)[refBaseMac],\
         (qrySeqStr)[qryBaseMac],\
         gapDiffS,\
         (scoreAry)[refBaseMac],\
         (dirAry)[refBaseMac],\
         macro_score_for_next_snp_L,\
         macroDelScoreL,\
         (alnSetPtr)\
      ); /*Get max score for the end of the row*/\
      \
      /*Prepare for the next round*/\
      macro_score_for_next_snp_L = macroGapColL;\
      macroGapColL += (alnSetPtr)->gapExtendC;\
      \
      macroDelScoreL =\
         macroGapColL + (alnSetPtr)->gapExtendC;\
   } /*Loop: score all query bases (rows)*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-04:
   ^  - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*Correct for being on the next row*/\
   macroGapColL -= (alnSetPtr)->gapExtendC;\
   \
   macroGapColL; /*Macro equviulant to a return*/\
}) /*scoreHirschFor*/

/*-------------------------------------------------------\
| Fun-04: scoreHirschRev
|  - Scores the input query and reference sequence
|    backwards
| Input:
|  - refSeqStr:
|    o reference sequence to align
|  - refStart:
|    o Position of the first base to align in the
|      reference sequence (is index 0)
|  - refLen:
|    o Number of bases to align in the reference sequence
|      (is index 1)
|  - qrySeqStr:
|    o query sequence to align
|  - qryStart:
|    o Position of the first base to align in the query
|      sequence (is index 0)
|  - qryLen:
|    o Number of bases to align in the query sequence
|      (is index 1)
|  - scoreAry:
|    o Array to fill with scores
|  - dirAry:
|    o array to use in finding directions
|  - alnSetPtr:
|    o Pointer to alnSet structure with the settingsq for
|      the alignment
| Output:
|  - Modifies:
|    o scoreAry to hold the last row of scores in a
|      Needleman Wunsch
|    o dirAry to hold the last row of directions in a
|      Needleman Wunsch
|  - Returns:
|    o The gap column score
\-------------------------------------------------------*/
#define scoreHirschRev(\
  refSeqStr,       /*Reference sequence*/\
  refStart,        /*index 0 starting ref base*/\
  refLen,          /*index 1 Length of ref*/\
  \
  qrySeqStr,       /*Query sequence*/\
  qryStart,        /*Index 0 Starting query base*/\
  qryLen,          /*index 1 length of query*/\
  \
  scoreAry,        /*Array of scores to fill*/\
  dirAry,          /*direction row, for gap extend*/\
  alnSetPtr        /*setttings to use*/\
)({/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: scoreHirschRev
   '  - Does a single round of scoring for a hirschberg
   '    alignment (backwards direction)
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-04 sec-03:
   '    - Score till on the last row
   '  o fun-04 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   short gapDiffS =\
      (alnSetPtr)->gapExtendC - (alnSetPtr)->gapOpenC;\
      /*Used to find the gap score (branchless method)*/\
   \
   ulong refBaseMac = 0;\
   ulong qryBaseMac = 0;\
   long macroGapColL = 0;\
   long macro_score_for_next_snp_L = 0;\
   long macroDelScoreL = 0;\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Set up the first row (indel row) of scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   (dirAry)[refStart + refLen - 1] = -1;\
   macroGapColL = 0;\
   (scoreAry)[refStart+refLen-1] = (alnSetPtr)->gapOpenC;\
   \
   for(\
      refBaseMac = (refStart) + (refLen) - 2;\
      refBaseMac > (refStart);\
      --refBaseMac\
   ){ /*Loop:Set the initial blank scores*/\
      (dirAry)[refBaseMac] = -1;\
      \
      (scoreAry)[refBaseMac] =\
           (scoreAry)[refBaseMac + 1]\
         + (alnSetPtr)->gapExtendC;\
   } /*Loop:Set the initial blank scores*/\
   \
   (dirAry)[refBaseMac] = -1;\
   (scoreAry)[refBaseMac] =\
      (scoreAry)[refStart + 1] + (alnSetPtr)->gapExtendC;\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Find the scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*Get the value I need to find the next snp's score*/\
   macro_score_for_next_snp_L = 0;\
   macroGapColL = (alnSetPtr)->gapOpenC;\
   \
   macroDelScoreL = macroGapColL + alnSetPtr->gapExtendC;\
   \
   for(\
      qryBaseMac = (qryStart) + (qryLen) - 1;\
      qryBaseMac >= (qryStart);\
      --qryBaseMac\
   ){ /*Loop: score all query bases (rows)*/\
      \
      for(\
         refBaseMac = (refStart) + (refLen) - 1;\
         refBaseMac > (refStart);\
         --refBaseMac\
      ){ /*Loop: Find the max scores for a single row*/\
         hirschScore(\
            (refSeqStr)[refBaseMac],\
            (qrySeqStr)[qryBaseMac],\
            gapDiffS,\
            (scoreAry)[refBaseMac],\
            (dirAry)[refBaseMac],\
            macro_score_for_next_snp_L,\
            macroDelScoreL,\
            (alnSetPtr)\
         ); /*Get max score for the current base pair*/\
      } /*Loop: Find the max scores for a single row*/\
      \
      hirschScoreRowEnd(\
         (refSeqStr)[refBaseMac],\
         (qrySeqStr)[qryBaseMac],\
         gapDiffS,\
         (scoreAry)[refBaseMac],\
         (dirAry)[refBaseMac],\
         macro_score_for_next_snp_L,\
         macroDelScoreL,\
         (alnSetPtr)\
      ); /*Get max score for the end of the row*/\
      \
      /*Prepare for the next round*/\
      macro_score_for_next_snp_L = macroGapColL;\
      macroGapColL += (alnSetPtr)->gapExtendC;\
      \
      macroDelScoreL =\
            macroGapColL + (alnSetPtr)->gapExtendC;\
   } /*Loop: score all query bases (rows)*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-04:
   ^  - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*Correct for being on the next row*/\
   macroGapColL -= (alnSetPtr)->gapExtendC;\
   \
   macroGapColL; /*Macro equviulant to a return*/\
}) /*scoreHirschRev*/

#endif
