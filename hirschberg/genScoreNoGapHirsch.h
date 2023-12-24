/*########################################################
# Name: generalScoreHirschNoGap
# Use:
#  - Holds functions for doing scoring in a scalar
#    Hirschberg alignment without a gap extension penalty
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
'  o fun-01 hirschScoreNoGap:
'    - Gets a single score for a single base pair in an
'      Hirschberg alignment without gap extensions penatly
'  o fun-02 hirschScoreRowEndNoGap:
'    - Finds the last score for the last base pair in a
'      row for a Hirschberg alignment without a gap
'      extension penalty
'  o fun-03 scoreHirschForNoGap:
'    - Does the scoring step whithout gap extensions
'      penalties for a hirschberg alignment
'      (forward direction)
'  o fun-04 scoreHirschRevNoGap:
'    - Does the scoring step whithout gap extensions
'      penalties for a hirschberg alignment
'      (reverse direction)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Includes and definitions
\-------------------------------------------------------*/

#ifndef GENERAL_SCORE_HIRSCH_NO_GAP_H 
#define GENERAL_SCORE_HIRSCH_NO_GAP_H 

#include "genHirsch.h"

/*-------------------------------------------------------\
| Fun-01: hirschScoreNoGap
|   - Gets a single score for a single base pair in an
|     Hirschberg alignment without gap extensions penatly
| Input:
|   - refBase:
|     o Reference base to get the score for (character)
|   - qryBase:
|     o Query base to get the score for (character)
|   - scoreOn:
|     o Score to updated and used to find insertions
|   - nextSnpScore:
|     o Has the score for finding an snp.
|     o This is updated to hold the next snps score
|   - delScore:
|     o Score for a deletion
|     o This is updated to hold the next deletion
|   - alnSetPtr:
|     o Pointer to an alnSet structer with settings to
|       use in the alignment
| Output:
|   - Modifies:
|     o scoreOn to hold the next maximum score
|     o nextSnpScore to hold the next score needed to
|       calculate the next snp
|     o delScore to hold the next deletion score
\-------------------------------------------------------*/
#define hirschScoreNoGap(\
   refBase,  /*Reference base*/\
   qryBase,  /*Query sequence*/\
   scoreOn,  /*Score on and holds the new maximum*/\
   nextSnpScore, /*Holds score used to get snp score*/\
   delScore,  /*Holds the deletion score*/\
   alnSetPtr  /*Holds the settings for the alignment*/\
){ /*hirschScoreNoGap*/\
   long macroSnpScoreL =\
        (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
   \
   /*Get the insertion score*/\
   long macroInsScoreL=(scoreOn) + (alnSetPtr)->gapOpenC;\
   \
   /*Get the score needed to find the next snp/match. This
   ` is overwiten when I find the max score.
   */\
   (nextSnpScore) = (scoreOn);\
   \
   alnMaxScore(\
      (scoreOn),\
      macroSnpScoreL,           /*snp/match score*/\
      macroInsScoreL,           /*insertion score*/\
      (delScore),               /*deletion score*/\
      (alnSetPtr)->bestDirC\
   );\
   \
   (delScore) = (scoreOn) + (alnSetPtr)->gapOpenC;\
} /*hirschScoreNoGap*/

/*-------------------------------------------------------\
| Fun-02: hirschScoreEnd
|   - Finds the last score for the last base pair in a
|     row for a Hirschberg alignment without a gap
|     extension penalty
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
|     o Pointer to an alnSet structer with settings to
|       use in the alignment
| Output:
|   - Modifies:
|     o scoreOn to hold the next maximum score
|     o dirOn to hold the next direction
\-------------------------------------------------------*/
#define hirschScoreRowEndNoGap(\
   refBase,  /*Reference base*/\
   qryBase,  /*Query sequence*/\
   scoreOn,  /*Holds new maxium score*/\
   nextSnpScore, /*Holds score used to get snp score*/\
   delScore,  /*Holds the deletion score*/\
   alnSetPtr  /*Holds the settings for the alignment*/\
){ /*hirschScoreRowEndNoGap*/\
   long macroSnpScoreL =\
        (nextSnpScore)\
      + getBaseScore((qryBase), (refBase), (alnSetPtr));\
   \
   /*Get the insertion score*/\
   long macroInsScoreL=(scoreOn) + (alnSetPtr)->gapOpenC;\
   \
   alnMaxScore(\
      (scoreOn),\
      macroSnpScoreL,           /*snp/match score*/\
      macroInsScoreL,           /*insertion score*/\
      (delScore),               /*deletion score*/\
      (alnSetPtr)->bestDirC\
   );\
} /*hirschScoreRowEndNoGap*/

/*-------------------------------------------------------\
| Fun-03: scoreHirschForNoGap
|  - Scores the input query and reference sequence for a
|    Hirscberg without a gap extension penalty
| Input:
|  - refSeqStr:
|    o C-string with reference sequence to align
|  - refStartUL:
|    o Position of the first base to align in the
|      reference sequence (is index 0)
|  - refLenUL:
|    o Number of bases to align in the reference sequence
|      (is index 1)
|  - qrySeqStr:
|    o C-string with query sequence to align
|  - qryStart:
|    o Position of the first base to align in the query
|      sequence (is index 0)
|  - qryLen:
|    o Number of bases to align in the query sequence
|      (is index 1)
|  - scoreRowPtrL:
|    o Array to fill with scores
|  - alnSetPtr:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      Needleman Wunsch
|  - Returns:
|    o The gap column score
\-------------------------------------------------------*/
#define scoreHirschForNoGap(\
  refSeqStr,      /*Reference sequence*/\
  refStart,       /*index 0 starting ref base*/\
  refLen,         /*index 1 Length of ref*/\
  \
  qrySeqStr,      /*Query sequence*/\
  qryStart,       /*Index 0 Starting query base*/\
  qryLen,         /*index 1 length of query*/\
  \
  scoreAry,       /*Array of scores to fill*/\
  alnSetPtr       /*setttings to use*/\
)({/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: scoreHirschForNoGap
   '  - Does the scoring for a hirschberg alignment
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
   \
   (scoreAry)[refStart] = (alnSetPtr)->gapOpenC;\
   \
   for(\
      refBaseMac = (refStart) + 1;\
      refBaseMac < (refStart) + (refLen);\
      ++refBaseMac\
   ){ /*Loop:Set the initial blank scores*/\
     (scoreAry)[refBaseMac] =\
        (scoreAry)[refBaseMac-1] + (alnSetPtr)->gapOpenC;\
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
   macroDelScoreL = macroGapColL + (alnSetPtr)->gapOpenC;\
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
         hirschScoreNoGap(\
            (refSeqStr)[refBaseMac],\
            (qrySeqStr)[qryBaseMac],\
            (scoreAry)[refBaseMac],\
            macro_score_for_next_snp_L,\
            macroDelScoreL,\
            (alnSetPtr)\
         ); /*Get max score for the current base pair*/\
      } /*Loop: Find the max scores for a single row*/\
      \
      hirschScoreRowEndNoGap(\
         (refSeqStr)[refBaseMac],\
         (qrySeqStr)[qryBaseMac],\
         (scoreAry)[refBaseMac],\
         macro_score_for_next_snp_L,\
         macroDelScoreL,\
         (alnSetPtr)\
      ); /*Get max score for the end of the row*/\
      \
      /*Prepare for the next round*/\
      macro_score_for_next_snp_L = macroGapColL;\
      \
      macroGapColL += (alnSetPtr)->gapOpenC;\
      macroDelScoreL=macroGapColL+(alnSetPtr)->gapOpenC;\
   } /*Loop: score all query bases (rows)*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-04:
   ^  - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*Correct for being on the next row*/\
   macroGapColL -= (alnSetPtr)->gapOpenC;\
   \
   macroGapColL; /*Macro equviulant to a return*/\
}) /*scoreHirschForNoGap*/

/*-------------------------------------------------------\
| Fun-04: scoreHirschRevNoGap
|  - Scores the input query and reference sequence
|    backwards. This function does not use the gap
|    extension penalty
| Input:
|  - refSeqStr:
|    o C-string with reference sequence to align
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
|  - scoreRowPtrL:
|    o Array to fill with scores
|  - alnSetPtr:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      Needleman Wunsch
|  - Returns:
|    o The gap column score
\-------------------------------------------------------*/
#define scoreHirschRevNoGap(\
  refSeqStr,      /*Reference sequence*/\
  refStart,       /*index 0 starting ref base*/\
  refLen,         /*index 1 Length of ref*/\
  \
  qrySeqStr,      /*Query sequence*/\
  qryStart,       /*Index 0 Starting query base*/\
  qryLen,         /*index 1 length of query*/\
  \
  scoreAry,       /*Array of scores to fill*/\
  alnSetPtr       /*setttings to use*/\
)({/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: scoreHirschRevNoGap
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
   macroGapColL = 0;\
   \
   (scoreAry)[refStart + refLen - 1] =\
      (alnSetPtr)->gapOpenC;\
   \
   for(\
      refBaseMac = (refStart) + (refLen) - 2;\
      refBaseMac > (refStart);\
      --refBaseMac\
   ){ /*Loop:Set the initial blank scores*/\
      (scoreAry)[refBaseMac] =\
        (scoreAry)[refBaseMac+1] + (alnSetPtr)->gapOpenC;\
   } /*Loop:Set the initial blank scores*/\
   \
   (scoreAry)[refBaseMac] =\
        (scoreAry)[refStart + 1] + (alnSetPtr)->gapOpenC;\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Find the scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*Get the value I need to find the next snp's score*/\
   macro_score_for_next_snp_L = 0;\
   macroGapColL = (alnSetPtr)->gapOpenC;\
   macroDelScoreL = macroGapColL + (alnSetPtr)->gapOpenC;\
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
         hirschScoreNoGap(\
            (refSeqStr)[refBaseMac],\
            (qrySeqStr)[qryBaseMac],\
            (scoreAry)[refBaseMac],\
            macro_score_for_next_snp_L,\
            macroDelScoreL,\
            (alnSetPtr)\
         ); /*Get the score for the base pair*/\
      } /*Loop: Find the max scores for a single row*/\
      \
      hirschScoreRowEndNoGap(\
         (refSeqStr)[refBaseMac],\
         (qrySeqStr)[qryBaseMac],\
         (scoreAry)[refBaseMac],\
         macro_score_for_next_snp_L,\
         macroDelScoreL,\
         (alnSetPtr)\
      ); /*Get max score for the end of the row*/\
      \
      /*Prepare for the next round*/\
      macro_score_for_next_snp_L = macroGapColL;\
      \
      macroGapColL += (alnSetPtr)->gapOpenC;\
      macroDelScoreL= macroGapColL+(alnSetPtr)->gapOpenC;\
   } /*Loop: score all query bases (rows)*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-04:
   ^  - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*Correct for being on the next row*/\
   macroGapColL -= (alnSetPtr)->gapOpenC;\
   \
   macroGapColL; /*Macro equviulant to a return*/\
}) /*scoreHirschRevNoGap*/

/*-------------------------------------------------------\
| Fun-05: positionSingleRefBaseNoGap
|  - Align a single base to a sequence
| Input:
|  - baseC:
|    o Base to align to a sequence
|  - baseIndexUL:
|    o Position of base in its seqequence (index 0)
|  - seqStr:
|    o Sequence to align the base to
|  - startOfSeqUL:
|    o Fist base to align in seqStr (index 0)
|  - lenSeqUL:
|    o Number of bases to align in seqStr (index 1)
|  - baseCAlnST:
|    o Holds the alignment of baseC's sequence
|    o if -DHIRSCHTWOBIT is two bit array, else char array
|  - seqAlnST:
|    o Holds the alignment of seqStr
|    o if -DHIRSCHTWOBIT is two bit array, else char array
|  - settings:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o baseCAlnST to hold the alignment for the single
|      base (baseC)
|    o seqAlnST to hold the alignment for the sequence
|      (seqStr)
\-------------------------------------------------------*/
static void positionSingleBaseNoGap(
  char baseC,        /*Single base to align to sequence*/
  ulong baseIndexUL, /*Index baseC is at*/
  char *seqStr,      /*Sequence to align baseC to*/
  ulong startOfSeqUL,/*1st base to align in seqStr*/
  ulong lenSeqUL,    /*number bases to align in seqStr*/
  char *baseAlnAryC, /*holds baseC alignment*/
  char *seqAlnAryC,  /*holds seqStr alignment*/
  struct alnSet *settings /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBaseNoGap
   '  - Align a single base to a sequence
   '  o fun-05 sec-01:
   '    - Variable declerations
   '  o fun-05 sec-02:
   '    - Find the reference bases position on the query
   '  o fun-05 sec-03:
   '    - Fill in insertions and reference base position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long insScoreL = 0;      /*The score of an insertion*/
   long matchScoreL = 0;    /*The score of a match*/

   ulong ulBase = 0;        /*for the for loop*/
   ulong snpIndexUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-02:
   ^  - Find the reference bases position on the query
   ^  o fun-05 sec-02 sub-01:
   ^    - Find the first scores for the loop
   ^  o fun-05 sec-02 sub-02:
   ^    - Find the remaing scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-05 Sec-02 Sub-01:
   *  - Find the first scores for the loop
   \*****************************************************/

   insScoreL = 0;

   /*****************************************************\
   * Fun-05 Sec-02 Sub-02:
   *  - Find the remaing scores
   \*****************************************************/
 
   for(
      ulBase = startOfSeqUL;
      ulBase < startOfSeqUL + lenSeqUL;
      ++ulBase
   ){ /*Loop: align baseC to seqStr*/
     matchScoreL =
         insScoreL
       + getBaseScore(seqStr[ulBase], baseC, settings);

     insScoreL = insScoreL + settings->gapOpenC;
     if(matchScoreL > insScoreL) snpIndexUL = ulBase;
   } /*Loop: align baseC to seqStr*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-03:
   ^  - Fill in the insertions and reference base position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   /* A series of deletions and insertions are prefered
   `  over matches and smps. In this case put the base
   `  at the first base. There is no good position
   */
   if(snpIndexUL == 0) snpIndexUL = startOfSeqUL;

   /*Add in the insertions at the start*/
   for(
      ulBase = startOfSeqUL;
      ulBase < snpIndexUL;
      ++ulBase
   ) seqAlnAryC[ulBase] = defGapFlag;
   
   /*Add in the position of the base*/
   baseAlnAryC[baseIndexUL] = defSnpFlag;
   seqAlnAryC[ulBase] = defSnpFlag;

   /*Finish adding in the insertions at the end*/
   for(
      ulBase = ulBase + 1;
      ulBase < startOfSeqUL + lenSeqUL;
      ++ulBase
   ) seqAlnAryC[ulBase] = defGapFlag;

   return;
} /*positionSingleRefBaseNoGap*/

#endif
