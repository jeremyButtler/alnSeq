/*#######################################################\
# Name: generalQryRefScan
#   - Holds functions to print and manage query refernce
#     scans 
# Libraries:
#   - "alnMatrixStruct.h" (print function) (No .c file)
#   - "genMath.h"                          (No .c file)
#   o "twoBitArrays.h"                     (No .c file)
#   o "dataTypeShortHand.h"                (No .c file)
# C Standard Libraries:
#   - <stdio.h>
\#######################################################*/


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  - Functions and variables for working with the query
'    reference scan Waterman alignments.
'  o header:
'    - definitions and includes
'  o fun-01 printAltWaterAlns:
'    - Prints out all saved alternatives alignments
'  o fun-02 getALnStart:
'    - Get the starting position, this is the old version
'  o set-01: Maximize score/direction (non-vector)
'    - Selects the max score and direction selected for
'      the max score.
'  o set-02: scanMaxScore
'    - This finds the max score based on the users
'      preferences for direction
'  o fun-03 waterScanMaxScore:
'    - Maximizes the score for a single base pair in an
'      waterman query reference scan
'  o fun-04 waterMaxEndRowScore:
'    - Maximizes the score for a single waterman alignment
'    - This for the end of a row, were you do not want the
'      next snp or deletion score
'  o fun-05 waterScanMaxScoreNoGap:
'    - Maximizes the score for a single base pair in an
'      waterman query reference scan. This function does
'      not apply a gap extension penalty
'  o fun-06 waterMaxEndRowScoreNoGap:
'    - Maximizes the score for a single waterman alignment
'      and does not apply a gap extension penalty for
'      deletions and insertions
'    - This for the end of a row, were you do not want the
'      next snp or deletion 
'  o set-03 scanIfKeepScore:
'    - Checks to see if the score should be kept in a
'      query reference scan. Each score is only recored
'      once.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*------------------------------------------------------\
| Header: 
|   - Definitions and includes
\------------------------------------------------------*/

#ifndef GENERAL_QRY_REF_SCAN_H
#define GENERAL_QRY_REF_SCAN_H

#include <stdio.h>

#include "alnMatrixStruct.h"
#include "genMath.h"

/*-------------------------------------------------------\
| Fun-01: printAltWaterAlns
|   - Prints out all saved alternatives alignments
| Input
|   - alnMtxPtr:
|     o Pointer to alnMatrix structure with alternative
|       alignments to print out
|   - minScore:
|     o Min score to keep an alternative alignment
|   - outFILE:
|     o File to write the alignment starting and ending
|       positions to
|   - refOffset:
|     o First base in reference alignment
|   - qryOffset:
|     o First base in query alignment
| Output:
|  - Prints
|    o Prints out the score, position of the first refence
|      base, last reference base, first query base, and
|      last query base to outFILE
\-------------------------------------------------------*/
#define pAltAlnScores(\
   alnMtxPtr,\
   minScore,\
   outFILE\
){\
  ulong ulScore = 0;\
  ulong lenAryUL = (alnMtxPtr)->lenArraysUL;\
  ulong lenUL = (alnMtxPtr)->lenRefUL;\
  \
  long *scoreAryL = (alnMtxPtr)->scoreAryL;\
  ulong *startAryUL = (alnMtxPtr)->startIndexAryUL;\
  ulong *endAryUL = (alnMtxPtr)->endIndexAryUL;\
  \
  ulong refStartUL = 0;\
  ulong refEndUL = 0;\
  \
  ulong qryStartUL = 0;\
  ulong qryEndUL = 0;\
  \
  /*Print out the alternative alignment header*/\
  fprintf(outFILE, "score\trefStart\trefEnd");\
  fprintf(outFILE, "\tqueryStart\tQueryEnd\n");\
  \
  for(ulScore = 0; ulScore < lenAryUL; ++ulScore)\
  { /*For all reference bases in the alignment*/\
    /*Check if score is worth printing out*/\
    if(scoreAryL[ulScore] < (minScore)) continue;\
    \
    indexToCoord(\
       lenUL,\
       startAryUL[ulScore],\
       refStartUL,\
       qryStartUL\
    );\
    \
    indexToCoord(\
       lenUL,\
       endAryUL[ulScore],\
       refEndUL,\
       qryEndUL\
    );\
    \
    fprintf(\
       outFILE,\
       "%li\t%lu\t%lu\t%lu\t%lu\n",\
       scoreAryL[ulScore],\
       refStartUL + (alnMtxPtr->refOffsetUL),\
       refEndUL + (alnMtxPtr->refOffsetUL),\
       qryStartUL + (alnMtxPtr->qryOffsetUL),\
       qryEndUL + (alnMtxPtr->qryOffsetUL)\
    ); /*Print out the coordinates*/\
  } /*For all reference bases in the alignment*/\
} /*pAltWaterAlns*/

/*-------------------------------------------------------\
| Fun-02: getALnStart
|   - Get the starting position, this is the old version
| Input:
|   - dir:
|     o Direction the current cell points to
|   - snpIndex:
|     o Start of the alignment if best dir was an match/snp
|   - insIndex:
|     o Start of the alignment if best dir was an insertion
|   - delIndex:
|     o Start of the alignment if best dir was an deletion
|   - newIndex:
|     o Index of dir in the alignment, it is used when
|       the direction was a stop.
| Output:
|  - Returns:
|    o The index of the input direction
\-------------------------------------------------------*/
/*The banchless version is slighty faster than using
' a switch statement
*/
#define getAlnStart(\
      dir,\
      snpIndex,\
      insIndex,\
      delIndex,\
      newIndex\
   )({\
      ulong retIndexUL = newIndex & -(dir == defMvStop);\
      retIndexUL += (snpIndex & -(dir == defMvSnp));\
      retIndexUL += (insIndex & -(dir == defMvIns));\
      retIndexUL += (delIndex & -(dir == defMvDel));\
      retIndexUL;\
   }) /*getAlnStart*/

/*-------------------------------------------------------\
| Set-01: Maximize score/direction (non-vector)
|  - Selects the max score and direction selected for the
|    max score.
|  o Set-01 a: insDelSnp (no braching)
|  o set-01 b: delInsSnp (no branching)
|  o set-01 c: insSnpDel (no branching)
|  o set-01 d: delSnpIns (no branching)
|  o set-01 e: snpInsDel (no branching)
|  o set-01 f: snpDelIns (no branching)
| Input:
|  - maxSc:
|    o This will hold the max score
|  - maxDir
|    o This will hold the direction of the max score
|  - maxPos:
|    o Index to store
|  - insSc:
|    o Score for having an insertion at this position
|  - snpSc
|    o Score for having an SNP/match at this position
|  - delSc:
|    o Score for having an deletion at this position
|  - snpPos:
|    o Index for and snp
|  - insPos:
|    o Index for a insertion
|  - delPos:
|    o Index for an deletion
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\-------------------------------------------------------*/

/*-----------------------------------------------------\
| Set-01 a: insDelSnp
\-----------------------------------------------------*/
#define scanInsDelSnp(\
   maxSc,\
   maxDir,\
   maxPos,\
   snpSc,\
   insSc,\
   delSc,\
   snpPos,\
   insPos,\
   delPos\
){\
   macroMax((maxSc), (insSc), (delSc));          /*5 Op*/\
   (maxDir) = (snpSc) > (maxSc);                 /*1 Op*/\
   macroIfMax(\
      (maxPos),\
      (insSc),\
      (delSc),\
      (insPos),\
      (delPos)\
   );\
   \
   macroIfMax(\
      (maxPos),\
      (maxSc),\
      (snpSc),\
      (maxPos),\
      (snpPos)\
   );\
   macroMax((maxSc), (maxSc), (snpSc));          /*5 Op*/\
   (maxDir) +=                                   /*4 Op*/\
      ( ( (insSc) == (maxSc) ) | (maxDir) ) + defMvDel;\
   /*Logic:
   ` maxDir = snp > max
   `   1 if snp was selected, else 0 (ins/del)
   ` defMvDel = 1
   ` (ins == max) | maxDir
   `   1 if insertion or snp was selected, else 0 (del)
   ` For an snp 1 += (0 | 1) + 1 = 3
   ` For an ins 0 += (1 | 0) + 1 = 2
   ` For an del 0 += (0 | 0) + 1 = 1
   */\
}

/*-----------------------------------------------------\
| Set-01 b: delInsSnp
\-----------------------------------------------------*/
#define scanDelInsSnp(\
   maxSc,\
   maxDir,\
   maxPos,\
   snpSc,\
   insSc,\
   delSc,\
   snpPos,\
   insPos,\
   delPos\
){\
   macroMax((maxSc), (delSc), (insSc));          /*5 Op*/\
   maxDir = (snpSc) > (maxSc);                   /*1 Op*/\
   macroIfMax(\
      (maxPos),\
      (delSc),\
      (insSc),\
      (delPos),\
      (insPos)\
   );\
   \
   macroIfMax(\
      (maxPos),\
      (maxSc),\
      (snpSc),\
      (maxPos),\
      (snpPos)\
   );\
   \
   macroMax((maxSc), (maxSc), (snpSc));          /*5 Op*/\
   (maxDir) += ( (delSc) != (maxSc) ) + defMvDel; /*3*/\
   /*Logic:
   ` maxDir = snp > max
   `   maxDir is 1 if snp was selected, else 0
   ` defMvDel = 1
   ` del != maxSc
   `   1 if a deletion (ins/snp) was selected
   ` For an snp 1 += (1) + 1 = 3
   ` For an ins 0 += (1) + 1 = 2
   ` For an del 0 += 0 + 1 = 1
   */\
}

/*-----------------------------------------------------\
| Set-01 c: insSnpDel
\-----------------------------------------------------*/
#define scanInsSnpDel(\
   maxSc,\
   maxDir,\
   maxPos,\
   snpSc,\
   insSc,\
   delSc,\
   snpPos,\
   insPos,\
   delPos\
){\
   macroMax(maxSc, (insSc), (snpSc));            /*5 Op*/\
   (maxDir) = (delSc) <= (maxSc);                /*1 Op*/\
     /*1 if kept score is not a deletion, else 0*/\
   macroIfMax(\
      (maxPos),\
      (insSc),\
      (snpSc),\
      (insPos),\
      (snpPos)\
   );\
   \
   macroIfMax(\
      (maxPos),\
      (maxSc),\
      (delSc),\
      (maxPos),\
      (delPos)\
   );\
   \
   macroMax((maxSc), (maxSc), (delSc));          /*5 Op*/\
   (maxDir) += (((snpSc) > (insSc)) & maxDir) + defMvDel;\
      /*4 Op*/\
   /*Logic:
   ` maxDir del <= max
   `   1 if an insertion or snp, else 0 (deletion)
   ` defMvDel = 1
   ` (snp > ins) & maxDir
   `   1 if an snp is selected, else 0 (ins/del)
   ` For an snp 1 += (1 & 1) + 1 = 3
   ` For an ins 1 += (0 & 1) + 1 = 2
   ` For an del 0 += (0/1 & 0) + 1 = 1
   */\
}

/*-----------------------------------------------------\
| Set-01 d: delSnpIns
\-----------------------------------------------------*/
#define scanDelSnpIns(\
   maxSc,\
   maxDir,\
   maxPos,\
   snpSc,\
   insSc,\
   delSc,\
   snpPos,\
   insPos,\
   delPos\
){\
   macroMax((maxSc), (snpSc), (insSc));          /*5 Op*/\
   (maxDir) = (maxSc) <= (delSc);                /*1 Op*/\
   macroIfMax(\
      (maxPos),\
      (snpSc),\
      (insSc),\
      (snpPos),\
      (insPos)\
   );\
   macroIfMax(\
      (maxPos),\
      (delSc),\
      (maxSc),\
      (delPos),\
      (maxPos)\
   );\
   macroMax((maxSc), (delSc), (maxSc));          /*5 Op*/\
   (maxDir) =                                    /*4 Op*/\
       defMvSnp\
     - ( ( (insSc) > (snpSc) ) | (maxDir) )\
     - (maxDir);\
   /*Logic:
   `  maxDir is 1 if a deletion is selected (max <= del)
   `  defMvSnp = 3
   `  (ins > snp) | (maxDir)
   `     1 if insertion or deletion selected, else 0
   `  -maxDir
   `     Is 1 if deletion was selected
   ` For an del I get 3 - (0/1 | 1) - 1 = 1
   ` For an ins I get 3 - (1 | 0) - 0 = 2
   ` For an snp I get 3 - (0 | 0) - 0 = 3
   */\
}

/*-----------------------------------------------------\
| Set-01 e: snpInsDel
\-----------------------------------------------------*/
#define scanSnpInsDel(\
   maxSc,\
   maxDir,\
   maxPos,\
   snpSc,\
   insSc,\
   delSc,\
   snpPos,\
   insPos,\
   delPos\
){\
   macroMax((maxSc), (snpSc), (insSc));          /*5 Op*/\
   (maxDir) = ((delSc) <= (maxSc));              /*1 Op*/\
     /*1 if deletion not choosen, else 0*/\
   macroIfMax(\
      (maxPos),\
      (snpSc),\
      (insSc),\
      (snpPos),\
      (insPos)\
   );\
   macroIfMax(\
      (maxPos),\
      (maxSc),\
      (delSc),\
      (maxPos),\
      (delPos)\
   );\
   macroMax((maxSc), (maxSc), (delSc));          /*5 Op*/\
   (maxDir) += ((snpSc) == (maxSc)) + defMvDel;  /*3 Op*/\
   /*Logic
   ` maxDir = delSc <= maxSc
   `   1 if an insertion or snp was selected, else 0 (del)
   ` snp == max
   `   1 if the snp was selected.
   ` For an snp 1 += (1) + 1 = 3
   ` For an ins 1 += (0) + 1 = 2
   ` For an del 0 += (0) + 1 = 1
   */\
}

/*-----------------------------------------------------\
| Set-01 f: snpDelIns
\-----------------------------------------------------*/
#define scanSnpDelIns(\
   maxSc,\
   maxDir,\
   maxPos,\
   snpSc,\
   insSc,\
   delSc,\
   snpPos,\
   insPos,\
   delPos\
){\
   macroMax((maxSc), (delSc), (insSc));          /*5 Op*/\
   maxDir = (snpSc) < (maxSc);                   /*1 Op*/\
   macroIfMax(\
      (maxPos),\
      (delSc),\
      (insSc),\
      (delPos),\
      (insPos)\
   );\
   macoIfMax(\
      (maxPos),\
      (maxSc),\
      (snpSc),\
      (maxPos),\
      (snpPos)\
   );\
   macroMax((maxSc), (snpSc), (maxSc));          /*5 Op*/\
   (maxDir) =                                    /*4 Op*/\
        defMvSnp /*Always 3*/\
      - (maxDir) \
      - ( ( (delSc) == (maxSc) ) & (maxDir) );\
     /*Logic
     ` maxDir = snp < max
     `   1 if deletion/insertion selected, else 0 (snp)
     ` defMvSnp is always 3
     ` (del == max) & maxDir
     `   1 if a deletin was selected, else 0 (snp/ins)
     ` For an snp (over ins/snp)   3 - 0 - (0/1 & 0) = 3
     ` For an ins (ignore if eqal) 3 - 1 - (0 & 0) = 2
     ` For an del (over ins)       3 - 1 - (1 & 1) = 1
     */\
}

/*--------------------------------------------------------\
| Set-02: scanMaxScore
|  - This finds the max score based on the users
|    preferences for direction
| o set-02 a: Prioirty ins/del/snp
| o set-02 b: Prioirty del/ins/snp
| o set-02 c: Prioirty ins/snp/del
| o set-02 d: Prioirty del/snp/ins
| o set-02 e: Prioirty snp/ins/del
| o set-02 f: Prioirty snp/del/ins
| o set-02 g: User selects
| Input:
|  - maScore:
|    o This will hold the maximum score
|  - maxDir:
|    o This will hold the direction of the max score
|  - maxPos:
|    o Index to store
|  - snp:
|    o Score for making an snp/match move
|  - ins:
|    o Score for making an insertion move
|  - del:
|    o Score for making an deletion move
|  - snpPos:
|    o Index for and snp
|  - insPos:
|    o Index for a insertion
|  - delPos:
|    o Index for an deletion
|  - pref:
|    o Direction combination that is prefered when scores
|      are equal. This can be hardcoded in.
| Output:
|  - Modifies
|    o maxScore to hold the maxium score
|    o maxDir to hold the direction of the max score
\--------------------------------------------------------*/

/*------------------------------------------------------\
| Set-02 a: InsDelSnp
\------------------------------------------------------*/
#ifdef INSDELSNP
   #define scanMaxScore(\
      maxScore,\
      maxDir,\
      maxPos,\
      snp,\
      ins,\
      del,\
      snpPos,\
      insPos,\
      delPos,\
      pref\
   ){\
      scanInsDelSnp(\
        (maxScore),\
        (maxDir),\
        (maxPos),\
        (snp),\
        (ins),\
        (del),\
        (snpPos),\
        (insPos),\
        (delPos)\
      );\
   }

/*------------------------------------------------------\
| Set-02 b: DelInsSnp
\------------------------------------------------------*/
#elif defined DELINSSNP
   #define scanMaxScore(\
      maxScore,\
      maxDir,\
      maxPos,\
      snp,\
      ins,\
      del,\
      snpPos,\
      insPos,\
      delPos,\
      pref\
   ){\
      scanDelInsSnp(\
        (maxScore),\
        (maxDir),\
        (maxPos),\
        (snp),\
        (ins),\
        (del),\
        (snpPos),\
        (insPos),\
        (delPos)\
      );\
   }

/*------------------------------------------------------\
| Set-02 c: InsSnpDel
\------------------------------------------------------*/
#elif defined INSSNPDEL 
   #define scanMaxScore(\
      maxScore,\
      maxDir,\
      maxPos,\
      snp,\
      ins,\
      del,\
      snpPos,\
      insPos,\
      delPos,\
      pref\
   ){\
      scanInsSnpDel(\
        (maxScore),\
        (maxDir),\
        (maxPos),\
        (snp),\
        (ins),\
        (del),\
        (snpPos),\
        (insPos),\
        (delPos)\
      );\
   }

/*------------------------------------------------------\
| Set-02 d: DelSnpIns
\------------------------------------------------------*/
#elif defined DELSNPINS
   #define scanMaxScore(\
      maxScore,\
      maxDir,\
      maxPos,\
      snp,\
      ins,\
      del,\
      snpPos,\
      insPos,\
      delPos,\
      pref\
   ){\
      scanDelSnpIns(\
        (maxScore),\
        (maxDir),\
        (maxPos),\
        (snp),\
        (ins),\
        (del),\
        (snpPos),\
        (insPos),\
        (delPos)\
      );\
   }

/*------------------------------------------------------\
| Set-02 e: SnpInsDel
\------------------------------------------------------*/
#elif defined SNPINSDEL
   #define scanMaxScore(\
      maxScore,\
      maxDir,\
      maxPos,\
      snp,\
      ins,\
      del,\
      snpPos,\
      insPos,\
      delPos,\
      pref\
   ){\
      scanSnpInsDel(\
        (maxScore),\
        (maxDir),\
        (maxPos),\
        (snp),\
        (ins),\
        (del),\
        (snpPos),\
        (insPos),\
        (delPos)\
      );\
   }

/*------------------------------------------------------\
| Set-02 f: SnpDelIns
\------------------------------------------------------*/
#elif defined SNPDELINS
   #define scanMaxScore(\
      maxScore,\
      maxDir,\
      maxPos,\
      snp,\
      ins,\
      del,\
      snpPos,\
      insPos,\
      delPos,\
      pref\
   ){\
      scanSnpDelIns(\
        (maxScore),\
        (maxDir),\
        (maxPos),\
        (snp),\
        (ins),\
        (del),\
        (snpPos),\
        (insPos),\
        (delPos)\
      );\
   }

/*------------------------------------------------------\
| Set-02 g: User selects
\------------------------------------------------------*/
#else
   #define scanMaxScore(\
      maxScore,\
      maxDir,\
      maxPos,\
      snp,\
      ins,\
      del,\
      snpPos,\
      insPos,\
      delPos,\
      pref\
   ){\
      switch(pref)\
      { /*Switch; get an snp/match priority*/\
         case defInsDelSnp:\
             scanInsDelSnp(\
               (maxScore),\
               (maxDir),\
               (maxPos),\
               (snp),\
               (ins),\
               (del),\
               (snpPos),\
               (insPos),\
               (delPos)\
             );\
            break;\
         case defDelInsSnp:\
             scanDelInsSnp(\
               (maxScore),\
               (maxDir),\
               (maxPos),\
               (snp),\
               (ins),\
               (del),\
               (snpPos),\
               (insPos),\
               (delPos)\
             );\
             break;\
         case defInsSnpDel:\
             scanInsSnpDel(\
               (maxScore),\
               (maxDir),\
               (maxPos),\
               (snp),\
               (ins),\
               (del),\
               (snpPos),\
               (insPos),\
               (delPos)\
             );\
             break;\
         case defDelSnpIns:\
             scanDelSnpIns(\
               (maxScore),\
               (maxDir),\
               (maxPos),\
               (snp),\
               (ins),\
               (del),\
               (snpPos),\
               (insPos),\
               (delPos)\
             );\
             break;\
         case defSnpInsDel:\
             scanSnpInsDel(\
               (maxScore),\
               (maxDir),\
               (maxPos),\
               (snp),\
               (ins),\
               (del),\
               (snpPos),\
               (insPos),\
               (delPos)\
             );\
             break;\
         case defSnpDelIns:\
             scanSnpDelIns(\
               (maxScore),\
               (maxDir),\
               (maxPos),\
               (snp),\
               (ins),\
               (del),\
               (snpPos),\
               (insPos),\
               (delPos)\
             );\
             break;\
      } /*Switch; get an snp/match priority*/\
   } /*charMaxScore*/
#endif /*For charMaxScore variations*/

/*-------------------------------------------------------\
| Fun-03: waterScanMaxScore
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
#define waterScanMaxScore(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,    /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
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
       + ((gapDiff) & (-((dirOn) != defMvSnp)))\
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
| Fun-04: waterMaxEndRowScore
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
#define waterScanMaxEndRowScore(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   gapDiff,    /*GapExtend - gapOpen*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
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
       + ((gapDiff) & (-((dirOn) != defMvSnp)))\
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

/*-------------------------------------------------------\
| Fun-05: waterScanMaxScoreNoGap
|   - Maximizes the score for a single base pair in an
|     waterman query reference scan. This function does
|     not apply a gap extension penalty
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
#define waterScanMaxScoreNoGap(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
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
    ulong keepDirUL = 0;\
    (nextSnpScore) = (scoreOn);/*Score to find next snp*/\
    \
    scanMaxScore(\
       (scoreOn),\
       (dirOn),\
       (index),\
       (macroSnpScoreL),\
       ((scoreOn) + (alnSetPtr)->gapOpenC),\
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
    (delScore) = (scoreOn) + (alnSetPtr)->gapOpenC;\
} /*waterScanMaxScoreNoGap*/

/*-------------------------------------------------------\
| Fun-06: waterMaxEndRowScoreNoGap
|   - Maximizes the score for a single waterman alignment
|     and does not apply a gap extension penalty for
|     deletions and insertions
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
#define waterScanMaxEndRowScoreNoGap(\
   refBase,     /*Reference Base*/\
   qryBase,     /*Query Base*/\
   scoreOn,     /*Score to update*/\
   dirOn,       /*Direction to update*/\
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
    ulong keepDirUL = 0;\
    \
    scanMaxScore(\
       (scoreOn),\
       (dirOn),\
       (index),\
       (macroSnpScoreL),\
       ((scoreOn) + (alnSetPtr)->gapOpenC),\
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

/*------------------------------------------------------\
| Set-03: scanIfKeepScore
|   - Checks to see if the score should be kept in a
|     query reference scan. Each score is only recored
|     once.
|   o set-03 a: scanIfKeepScoreRef
|     - Prioritizes the reference first and then the
|       query if the reference did not keep the score
|   o set-03 b: scanIfKeepScoreQry
|     - Prioritizes the query first and then the
|       reference if the query did not keep the score
| Input:
|   - score:
|     o Score to compare query and reference high scores
|       to
|   - start:
|     o Starting index to keep if score was better
|   - end:
|     o ending index to keep fi the score was better
|   - refScore:
|     o Current best score for the reference base
|   - refStart:
|     o Starting index of the reference alignment
|   - refEnd:
|     o Ending position of the reference alignment
|   - qryScore:
|     o Current best score for the query base
|   - qryStart:
|     o Starting index of the query alignment
|   - qryEnd:
|     o Ending position of the query alignment
| Output:
|   - Modifies:
|     o refScore to hold score if score > refScore
|     o refStart to hold start if score > refScore
|     o refEnd to hold end if score > refScore
|     o qryScore to hold score if score > qryScore
|     o qryStart to hold start if score > qryScore
|     o qryEnd to hold end if score > qryScore
\------------------------------------------------------*/

/*------------------------------------------------------\
| Set-03 a: scanIfKeepScoreRef
|   - Prioritizes the reference first and then the
|     query if the reference did not keep the score
\------------------------------------------------------*/
#define scanIfKeepScoreRef(\
   score,\
   dir,\
   start,\
   end,\
   refScore,\
   refStart,\
   refEnd,\
   qryScore,\
   qryStart,\
   qryEnd\
){\
   ulong oldRefBl;\
   ulong oldQryBl;\
   \
   ulong qryBl = ((dir) == defMvSnp);\
   ulong refBl = -(qryBl & ((refScore) < (score)));\
   qryBl = -((!refBl) & (qryBl & ((qryScore) < (score))));\
   \
   oldQryBl = ~qryBl;\
   oldRefBl = ~refBl;\
   \
   (refScore)=((refScore) & oldRefBl) + ((score) &refBl);\
   (refStart)=((refStart) & oldRefBl) + ((start) &refBl);\
   (refEnd) = ((refEnd) & oldRefBl) + ((end) & refBl);\
   \
   (qryScore)=((qryScore) & oldQryBl) + ((score) &qryBl);\
   (qryStart)=((qryStart) & oldQryBl) + ((start) &qryBl);\
   (qryEnd) = ((qryEnd) & oldQryBl) + ((end) & qryBl);\
} /*scanIfKeepScoreQry*/

/*------------------------------------------------------\
| Set-03 b: scanIfKeepScoreQry
|   - Prioritizes the query first and then the
|     query if the query did not keep the score
\------------------------------------------------------*/
#define scanIfKeepScoreQry(\
   score,\
   dir,\
   start,\
   end,\
   refScore,\
   refStart,\
   refEnd,\
   qryScore,\
   qryStart,\
   qryEnd\
){\
   ulong oldRefBl;\
   ulong oldQryBl;\
   \
   ulong refBl = (dir) == defMvSnp;\
   ulong qryBl = -(refBl & ((qryScore) < (score)));\
   refBl = -((!qryBl) & (refBl & ((refScore) < (score))));\
   \
   oldQryBl = ~qryBl;\
   oldRefBl = ~refBl;\
   \
   (qryScore)=((qryScore) & oldQryBl) + ((score) &qryBl);\
   (qryStart)=((qryStart) & oldQryBl) + ((start) &qryBl);\
   (qryEnd) = ((qryEnd) & oldQryBl) + ((end) & qryBl);\
   \
   (refScore)=((refScore) & oldRefBl) + ((score) &refBl);\
   (refStart)=((refStart) & oldRefBl) + ((start) &refBl);\
   (refEnd) = ((refEnd) & oldRefBl) + ((end) & refBl);\
} /*scanIfKeepScoreQry*/

#endif
