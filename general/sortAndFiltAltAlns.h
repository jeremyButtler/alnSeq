/*########################################################
# Name: scoreST
# Use:
#  - Has functions to sort and filter the arrays in the
#    alnMatrix or alnMatrixTwoBit structures
# Libraries:
#  - "alnMatrixStruct.h"    (No .c file)
#  o dataTypeShortHand.h"   (No .c file)
# C Standard Libraries:
#  o <stdlib.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - Includes and definitions
'  o fun-01 swap:
'    - Swaps values in two score structures
'    - Macro (in scoresST.h)
'  o fun-02: alnMatrixSortScores
'     - Sorts an array of scores structs by score using
'       shell short.  Order is greatest to least.
'  o fun-03: sortQryScores
'     - Sorts the three arrays in an alnMatrix structure
'       by the query starting position (lowest first) and
'       then by scores (highest first).
'  o fun-04: sortRefScores
'     - Sorts the three arrays in an alnMatrix structure
'       by the reference starting position (lowest first)
'       and then by scores (highest first).
'  o fun-05: sortRefQryScores
'     - Sorts the three arrays in an alnMatrix structure
'       by the query starting position (lowest first),
'       then the reference starting position (lowest
'       first), and then by the score (highest first).
'  o fun-06: altScoreFlitQry
'     - Filters by removing alternative alignments that
'       have an overlap in their position on the query.
'       The overlap with the highest score is kept.
'  o fun-07: altScoreFlitRef
'     - Filters by removing alternative alignments that
'       have an overlap in their position on the
'       reference The overlap with the highest score is
'       kept.
'  o fun-08: altScoreFlitRefQry
'     - Filters by removing alternative alignments that
'       have an overlap in their position on the reference
'       or query. The overlap with the highest score is
'       kept.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Includes and definitions
\-------------------------------------------------------*/

#ifndef SORT_AND_FILT_ALT_ALNS_H
#define SORT_AND_FILT_ALT_ALNS_H

#include "alnMatrixStruct.h"

/*-------------------------------------------------------\
| Fun-01: swapElm
|  - Swaps two (non-pointer) elements around
| Input:
|  - elmOne
|    o The first element to swap
|  - elmTwo
|    o The second element to swap
| Output:
|  - Modifies:
|    o elmOne position to hold elmTwo
|    o elmTwo position to hold elmOne
\-------------------------------------------------------*/
#define swap(elmOne, elmTwo){ \
  ulong macroSwapUL = 0;\
  \
  macroSwapUL = (elmOne);\
  (elmOne) = (elmTwo);\
  (elmTwo) = macroSwapUL;\
} /*swapScoreSTs*/

/*-------------------------------------------------------\
| Fun-02: alnMatrixSortScores
|  - Sorts an array of scores structs by score using
|    shell short.  Order is greatest to least.
| Input:
|  - mtrxST:
|    o Pointer to alnMatrix structure with scores to sort
|      and the indexes to keep track of.
|  - firstElmUL:
|    o First element to start sorting at
|  - endElmUL:
|    o Last element to sort
| Output:
|  - Modifies
|    o scoreAryL to be sorted by scores. This also makes
|      sure that startIndexAryUL and endIndexAryUL stay
|      in sync with scoreAryL.
\-------------------------------------------------------*/
#define alnMatrixSortScores(mtrxSTPtr, firstElm, endElm){\
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: alnMatrixSortScores
   '  - Sorts an array of scores structs by score using
   '    shell short.  Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-02 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
  \
  /*Number of elements to sort*/\
  ulong numElmUL = (endElm) - (firstElm);\
  \
  /*Number of sorting rounds*/\
  ulong subUL = 0;\
  ulong nextElmUL = 0;\
  ulong lastElmUL = 0;\
  ulong elmOnUL = 0;\
  \
  /*Get arrays to sort from the matrix (for sanity)*/\
  long *scoreAry = (mtrxSTPtr)->scoreAryL;\
  ulong *startAry = (mtrxSTPtr)->startIndexAryUL;\
  ulong *endAry = (mtrxSTPtr)->endIndexAryUL;\
  \
  /*Variables to incurment loops*/\
  ulong ulIndex = 0;\
  ulong ulElm = 0;\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\\
  ^ Fun-02 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] +1*/\
  subUL = 1; /*Initialzie first array*/\
  while(subUL < numElmUL) subUL = (3 * subUL) + 1;\
  subUL = (subUL - 1) / 3; /*Account for overshooting*/\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-02 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  while(subUL > 0)\
  { /*loop trhough all sub arrays sort the subarrays*/\
    for(ulIndex = 0; ulIndex <= subUL; ++ulIndex)\
    { /*For each element in the subarray*/\
      for(ulElm = ulIndex;\
          ulElm + subUL <= endElm;\
          ulElm += subUL\
      ){ /*Loop; swap each nth element of the subarray*/\
        nextElmUL = ulElm + subUL;\
        \
        if(scoreAry[ulElm] > scoreAry[nextElmUL])\
        { /*If I need to swap an element*/\
          swap(scoreAry[ulElm], scoreAry[nextElmUL]);\
          swap(startAry[ulElm], startAry[nextElmUL]);\
          swap(endAry[ulElm], endAry[nextElmUL]);\
          \
          lastElmUL = ulElm;\
          elmOnUL = ulElm;\
          \
          while(lastElmUL >= subUL)\
          { /*loop; move swapped element back*/\
            lastElmUL -= subUL;\
            \
            if(scoreAry[elmOnUL] > scoreAry[lastElmUL])\
               break; /*Positioned the element*/\
            \
            swap(scoreAry[elmOnUL], scoreAry[lastElmUL]);\
            swap(startAry[elmOnUL], startAry[lastElmUL]);\
            swap(endAry[elmOnUL], endAry[lastElmUL]);\
            \
            elmOnUL = lastElmUL;\
          } /*loop; move swapped element back*/\
        } /*If I need to swap elements*/\
      } /*Loop; swap each nth element of the subarray*/\
    } /*For each element in the subarray*/\
    \
    subUL = (subUL - 1) / 3; /*Move to the next round*/\
  } /*loop through all sub arrays to sort the subarrays*/\
} /*alnMatrixSortScores*/

/*-------------------------------------------------------\
| Fun-03: sortQryScores
|  - Sorts the three arrays in an alnMatrix structure by
|    the query starting position (lowest first) and then
|    scores (highest first).
| Input:
|  - mtrxST:
|    o Pointer to alnMatrix structure with scores to sort
|      and the indexes to keep track of.
|  - firstElmUL:
|    o First element to start sorting at
|  - endElmUL:
|    o Last element to sort
| Output:
|  - Modifies
|    o scoreAryL to be sorted by scores. This also makes
|      sure that startIndexAryUL and endIndexAryUL stay
|      in sync with scoreAryL.
\-------------------------------------------------------*/
#define sortQryScores(mtrxSTPtr, firstElm, endElm){\
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: srotQryScores
   '  - Sorts an array of scores structs by score using
   '    shell short.  Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-03 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
  \
  /*Number of elements to sort*/\
  ulong macroLenAryUL = (endElm) - (firstElm);\
  \
  /*Number of sorting rounds*/\
  ulong subUL = 0;\
  ulong nextElmUL = 0;\
  ulong lastElm = 0;\
  ulong elmOn = 0;\
  \
  /*Get arrays to sort from the matrix (for sanity)*/\
  long *scoreAry = (mtrxSTPtr)->scoreAryL;\
  ulong *startAry = (mtrxSTPtr)->startIndexAryUL;\
  ulong *endAry = (mtrxSTPtr)->endIndexAryUL;\
  \
  long scoreL = 0;\
  long nextScoreL = 0;\
  long lastScoreL = 0;\
  \
  ulong refLen = (mtrxSTPtr)->lenRefUL;\
  ulong qryPos = 0;\
  ulong nextQry = 0;\
  ulong lastQry = 0;\
  \
  /*Variables to incurment loops*/\
  ulong ulIndex = 0;\
  ulong ulElm = 0;\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\\
  ^ Fun-03 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] +1*/\
  subUL = 1; /*Initialzie first array*/\
  while(subUL < macroLenAryUL) subUL = (3 * subUL) + 1;\
  subUL = (subUL - 1) / 3; /*Account for overshooting*/\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  while(subUL > 0)\
  { /*loop trhough all sub arrays sort the subarrays*/\
    for(ulIndex = 0; ulIndex <= subUL; ++ulIndex)\
    { /*For each element in the subarray*/\
      for(ulElm = ulIndex;\
          ulElm + subUL <= endElm;\
          ulElm += subUL\
      ){ /*Loop; swap each nth element of the subarray*/\
        nextElmUL = ulElm + subUL;\
        \
        qryPos = indexToQry(refLen, startAry[ulElm]);\
        nextQry = indexToQry(refLen,startAry[nextElmUL]);\
        \
        scoreL = scoreAry[ulElm];\
        nextScoreL = scoreAry[nextElmUL];\
        \
        if(   qryPos < nextQry\
           || (qryPos == nextQry && scoreL > nextScoreL)\
        ){ /*If I need to swap an element*/\
          swap(scoreAry[ulElm], scoreAry[nextElmUL]);\
          swap(startAry[ulElm], startAry[nextElmUL]);\
          swap(endAry[ulElm], endAry[nextElmUL]);\
          \
          lastElm = ulElm;\
          elmOn = ulElm;\
          \
          while(lastElm >= subUL)\
          { /*loop; move swapped element back*/\
            lastElm -= subUL;\
            \
            qryPos = indexToQry(refLen, startAry[elmOn]);\
            lastQry=indexToQry(refLen,startAry[lastElm]);\
            \
            scoreL = scoreAry[elmOn];\
            lastScoreL = scoreAry[lastElm];\
            \
            if(   qryPos > lastQry\
               || (qryPos==lastQry && scoreL<=lastScoreL)\
            ) break; /*element in position*/\
            \
            swap(scoreAry[elmOn], scoreAry[lastElm]);\
            swap(startAry[elmOn], startAry[lastElm]);\
            swap(endAry[elmOn], endAry[lastElm]);\
            \
            elmOn = lastElm;\
          } /*loop; move swapped element back*/\
        } /*If I need to swap elements*/\
      } /*Loop; swap each nth element of the subarray*/\
    } /*For each element in the subarray*/\
    \
    subUL = (subUL - 1) / 3; /*Move to the next round*/\
  } /*loop through all sub arrays to sort the subarrays*/\
} /*sortQryScores*/

/*-------------------------------------------------------\
| Fun-04: sortRefScores
|  - Sorts the three arrays in an alnMatrix structure by
|    the reference starting position (lowest first) and
|    then scores (highest first).
| Input:
|  - mtrxST:
|    o Pointer to alnMatrix structure with scores to sort
|      and the indexes to keep track of.
|  - firstElmUL:
|    o First element to start sorting at
|  - endElmUL:
|    o Last element to sort
| Output:
|  - Modifies
|    o scoreAryL to be sorted by scores. This also makes
|      sure that startIndexAryUL and endIndexAryUL stay
|      in sync with scoreAryL.
\-------------------------------------------------------*/
#define sortRefScores(mtrxSTPtr, firstElm, endElm){\
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: sortRefScores
   '  - Sorts an array of scores structs by score using
   '    shell short.  Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-04 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
  \
  /*Number of elements to sort*/\
  ulong macroLenAryUL = (endElm) - (firstElm);\
  \
  /*Number of sorting rounds*/\
  ulong subUL = 0;\
  ulong nextElmUL = 0;\
  ulong lastElm = 0;\
  ulong elmOn = 0;\
  \
  /*Get arrays to sort from the matrix (for sanity)*/\
  long *scoreAry = (mtrxSTPtr)->scoreAryL;\
  ulong *startAry = (mtrxSTPtr)->startIndexAryUL;\
  ulong *endAry = (mtrxSTPtr)->endIndexAryUL;\
  \
  long scoreL = 0;\
  long nextScoreL = 0;\
  long lastScoreL = 0;\
  \
  ulong refLen = (mtrxSTPtr)->lenRefUL;\
  ulong refPos = 0;\
  ulong nextRef = 0;\
  ulong lastRef = 0;\
  \
  /*Variables to incurment loops*/\
  ulong ulIndex = 0;\
  ulong ulElm = 0;\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\\
  ^ Fun-04 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] +1*/\
  subUL = 1; /*Initialzie first array*/\
  while(subUL < macroLenAryUL) subUL = (3 * subUL) + 1;\
  subUL = (subUL - 1) / 3; /*Account for overshooting*/\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  while(subUL > 0)\
  { /*loop trhough all sub arrays sort the subarrays*/\
    for(ulIndex = 0; ulIndex <= subUL; ++ulIndex)\
    { /*For each element in the subarray*/\
      for(ulElm = ulIndex;\
          ulElm + subUL <= endElm;\
          ulElm += subUL\
      ){ /*Loop; swap each nth element of the subarray*/\
        nextElmUL = ulElm + subUL;\
        \
        refPos = indexToRef(refLen, startAry[ulElm]);\
        nextRef = indexToRef(refLen,startAry[nextElmUL]);\
        \
        scoreL = scoreAry[ulElm];\
        nextScoreL = scoreAry[nextElmUL];\
        \
        if(   refPos < nextRef\
           || (refPos == nextRef && scoreL > nextScoreL)\
        ){ /*If I need to swap an element*/\
          swap(scoreAry[ulElm], scoreAry[nextElmUL]);\
          swap(startAry[ulElm], startAry[nextElmUL]);\
          swap(endAry[ulElm], endAry[nextElmUL]);\
          \
          lastElm = ulElm;\
          elmOn = ulElm;\
          \
          while(lastElm >= subUL)\
          { /*loop; move swapped element back*/\
            lastElm -= subUL;\
            \
            refPos = indexToRef(refLen, startAry[elmOn]);\
            lastRef=indexToRef(refLen,startAry[lastElm]);\
            \
            scoreL = scoreAry[elmOn];\
            lastScoreL = scoreAry[lastElm];\
            \
            if(   refPos > lastRef\
               || (refPos==lastRef && scoreL<=lastScoreL)\
            ) break; /*element in position*/\
            \
            swap(scoreAry[elmOn], scoreAry[lastElm]);\
            swap(startAry[elmOn], startAry[lastElm]);\
            swap(endAry[elmOn], endAry[lastElm]);\
            \
            elmOn = lastElm;\
          } /*loop; move swapped element back*/\
        } /*If I need to swap elements*/\
      } /*Loop; swap each nth element of the subarray*/\
    } /*For each element in the subarray*/\
    \
    subUL = (subUL - 1) / 3; /*Move to the next round*/\
  } /*loop through all sub arrays to sort the subarrays*/\
} /*sortRefScores*/

/*-------------------------------------------------------\
| Fun-05: sortRefQryScores
|  - Sorts the three arrays in an alnMatrix structure by
|    the query starting position (lowest first), then the
|    reference starting position (lowest first), and then
|    the score (highest first).
| Input:
|  - mtrxST:
|    o Pointer to alnMatrix structure with scores to sort
|      and the indexes to keep track of.
|  - firstElmUL:
|    o First element to start sorting at
|  - endElmUL:
|    o Last element to sort
| Output:
|  - Modifies
|    o scoreAryL to be sorted by qry starting position,
|      reference starting position, and then scores. This
|      also makes sure that startIndexAryUL and
|      endIndexAryUL stay in sync with scoreAryL.
\-------------------------------------------------------*/
#define sortRefQryScores(mtrxSTPtr, firstElm, endElm){\
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: sortRefQryScores
   '  - Sorts an array of scores structs by score using
   '    shell short.  Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-05 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
  \
  /*Number of elements to sort*/\
  ulong macroNumElmUL = (endElm) - (firstElm);\
  \
  /*Number of sorting rounds*/\
  ulong subUL = 0;\
  ulong nextElmUL = 0;\
  ulong lastElmUL = 0;\
  ulong elmOnUL = 0;\
  \
  /*Get arrays to sort from the matrix (for sanity)*/\
  long *scoreAry = (mtrxSTPtr)->scoreAryL;\
  ulong *startAry = (mtrxSTPtr)->startIndexAryUL;\
  ulong *endAry = (mtrxSTPtr)->endIndexAryUL;\
  \
  ulong indexUL = 0;\
  ulong refLen = (mtrxSTPtr)->lenRefUL;\
  \
  long scoreL = 0;\
  long nextScoreL = 0;\
  long lastScoreL = 0;\
  \
  ulong qryPos = 0;\
  ulong nextQry = 0;\
  ulong lastQry = 0;\
  \
  ulong refPos = 0;\
  ulong nextRef = 0;\
  ulong lastRef = 0;\
  \
  /*Variables to incurment loops*/\
  ulong ulIndex = 0;\
  ulong ulElm = 0;\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\\
  ^ Fun-05 Sec-02:
  ^  - Find the max search value
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  /*macroNumElmUL = (endElm) - (firstElm);*/\
  \
  /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] +1*/\
  subUL = 1; /*Initialzie first array*/\
  while(subUL < macroNumElmUL) subUL = (3 * subUL) + 1;\
  subUL = (subUL - 1) / 3; /*Account for overshooting*/\
  \
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-05 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
  \
  while(subUL > 0)\
  { /*loop trhough all sub arrays sort the subarrays*/\
    for(ulIndex = 0; ulIndex <= subUL; ++ulIndex)\
    { /*For each element in the subarray*/\
      for(ulElm = ulIndex;\
          ulElm + subUL <= endElm;\
          ulElm += subUL\
      ){ /*Loop; swap each nth element of the subarray*/\
        nextElmUL = ulElm + subUL;\
        \
        indexUL = startAry[ulElm];\
        indexToCoord(refLen, indexUL, refPos, qryPos);\
        \
        indexUL = startAry[nextElmUL];\
        indexToCoord(refLen, indexUL, nextRef, nextQry);\
        \
        scoreL = scoreAry[ulElm];\
        nextScoreL = scoreAry[nextElmUL];\
        \
        if(   qryPos < nextQry\
           || refPos < nextRef\
           || (qryPos == nextQry && scoreL > nextScoreL)\
           || (refPos == nextRef && scoreL > nextScoreL)\
        ){ /*If I need to swap an element*/\
          swap(scoreAry[ulElm], scoreAry[nextElmUL]);\
          swap(startAry[ulElm], startAry[nextElmUL]);\
          swap(endAry[ulElm], endAry[nextElmUL]);\
          \
          lastElmUL = ulElm;\
          elmOnUL = ulElm;\
          \
          while(lastElmUL >= subUL)\
          { /*loop; move swapped element back*/\
            lastElmUL -= subUL;\
            \
            indexUL = startAry[elmOnUL];\
            indexToCoord(refLen,indexUL,refPos,qryPos);\
            \
            indexUL = startAry[lastElmUL];\
            indexToCoord(refLen,indexUL,lastRef,lastQry);\
            \
            scoreL = scoreAry[elmOnUL];\
            lastScoreL = scoreAry[lastElmUL];\
            \
            if(   qryPos > lastQry\
               || refPos > lastRef\
               || (qryPos==lastQry && scoreL<=lastScoreL)\
               || (refPos==lastRef && scoreL<=lastScoreL)\
            ) break; /*element in position*/\
            \
            swap(scoreAry[elmOnUL], scoreAry[lastElmUL]);\
            swap(startAry[elmOnUL], startAry[lastElmUL]);\
            swap(endAry[ulElm], endAry[nextElmUL]);\
            \
            elmOnUL = lastElmUL;\
          } /*loop; move swapped element back*/\
        } /*If I need to swap elements*/\
      } /*Loop; swap each nth element of the subarray*/\
    } /*For each element in the subarray*/\
    \
    subUL = (subUL - 1) / 3; /*Move to the next round*/\
  } /*loop through all sub arrays to sort the subarrays*/\
} /*sortRefQryScores*/

/*-------------------------------------------------------\
| Fun-06: altScoreFlitQry
|  - Filters by removing alternative alignments that
|    have an overlap in their position on the query.
|    The overlap with the highest score is kept.
| Input:
|  - mtrxSTPtr:
|    o Pointer to an alnMatrix or alnMatrixTwoBit
|      structure with scores to filter.
| Output:
|  - Modifies:
|    o scoreAryL, startIndexAryUL, and endIndexAryUL to
|      have overlapping scores removed
\-------------------------------------------------------*/
#define altScoreFiltQry(mtrxSTPtr){\
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: altScoreFlitQry
   '  - Filters by removing alternative alignments that
   '    have an overlap in their position on the query.
   '    The overlap with the highest score is kept.
   '  o fun-06 sec-01:
   '    - Variable declerations
   '  o fun-06 sec-02:
   '    - Remove overlapping alignments
   '  o fun-06 sec-03:
   '    - Blank the empty/discared alignments
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   ulong numElmUL = (mtrxSTPtr)->lenArraysUL;\
   ulong lenRefUL = (mtrxSTPtr)->lenRefUL;\
   \
   long *scorePtrL = (mtrxSTPtr)->scoreAryL;\
   long *lastScorePL = scorePtrL + numElmUL - 1;\
   ulong *startPtrUL = (mtrxSTPtr)->startIndexAryUL;\
   ulong *endPtrUL = (mtrxSTPtr)->endIndexAryUL;\
   \
   long *scoreOnPL = (mtrxSTPtr)->scoreAryL;\
   ulong *startOnPUL = (mtrxSTPtr)->startIndexAryUL;\
   ulong *endOnPUL = (mtrxSTPtr)->endIndexAryUL;\
   \
   ulong nextQryStartUL = 0;\
   ulong qryEndUL = 0;\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-02:
   ^  - Remove overlapping alignments
   ^  o fun-06 sec-02 sub-01:
   ^    - Loop setup + deal with blank or discared scores
   ^  o fun-06 sec-02 sub-02:
   ^    - Check if have an overlap
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*****************************************************\
   * Fun-06 Sec-02 Sub-01:
   *   - Loop setup + deal with blank or discared scores
   \****************************************************/\
   \
   sortQryScores((mtrxSTPtr), 0, numElmUL - 1);\
   \
   while(scoreOnPL <= lastScorePL)\
   { /*Loop: Check all scores to see if overlap*/\
      if(*scoreOnPL == 0)\
      { /*If: This had a score of 0 (was blank)*/\
         ++scoreOnPL;\
         ++startOnPUL;\
         ++endOnPUL;\
         continue;\
      } /*If: This had a score of 0 (was blank)*/\
      \
      if(*scorePtrL == 0)\
      { /*If: My comparision score is a blank*/\
         *scorePtrL = *scoreOnPL;\
         *startPtrUL = *startOnPUL;\
         *endPtrUL = *endOnPUL;\
         \
         ++scoreOnPL;\
         ++startOnPUL;\
         ++endOnPUL;\
         continue;\
      } /*If: My comparision score is a blank*/\
      \
      /**************************************************\
      * Fun-06 Sec-02 Sub-02:
      *   - Check if have an overlap
      \*************************************************/\
      \
      qryEndUL = indexToQry(lenRefUL, *endPtrUL);\
      nextQryStartUL = indexToQry(lenRefUL, *startOnPUL);\
      \
      if(nextQryStartUL <= qryEndUL)\
      { /*If: I have an overlap*/\
         if(*scorePtrL < *scoreOnPL)\
         { /*If: the new score is better*/\
            *scorePtrL = *scoreOnPL;\
            *startPtrUL = *startOnPUL;\
            *endPtrUL = *endOnPUL;\
         } /*If: the new score is better*/\
      } /*If: I have an overlap*/\
      \
      else if(scorePtrL != scoreOnPL)\
      { /*If: I am overwriting blank/overlapping score*/\
         ++scorePtrL;\
         ++startPtrUL;\
         ++endPtrUL;\
         *scorePtrL = *scoreOnPL;\
         *startPtrUL = *startOnPUL;\
         *endPtrUL = *endOnPUL;\
      } /*If: I am overwriting blank/overlapping score*/\
      \
      else\
      { /*Else: I need to move to the next score*/\
         ++scorePtrL;\
         ++startPtrUL;\
         ++endPtrUL;\
         /*case scorePtrL == scoreOnPL (same address)*/\
      } /*Else: I need to move to the next score*/\
      \
      ++scoreOnPL;\
      ++startOnPUL;\
      ++endOnPUL;\
   } /*Loop: Check all scores to see if overlap*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-03:
   ^  - Blank the empty/discared alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   if(scorePtrL < scoreOnPL)\
   { /*If: overwriting scores, get off last kept score*/\
      ++scorePtrL;\
      ++startPtrUL;\
      ++endPtrUL;\
   } /*If: overwriting scores, get off last kept score*/\
   \
   while(scorePtrL < scoreOnPL)\
   { /*Loop: blank ending scores*/\
      *scorePtrL = 0;\
      ++scorePtrL;\
   } /*Loop: blank ending scores*/\
} /*altScoreFiltQry*/

/*-------------------------------------------------------\
| Fun-07: altScoreFlitRef
|  - Filters by removing alternative alignments that
|    have an overlap in their position on the reference.
|    The overlap with the highest score is kept.
| Input:
|  - mtrxSTPtr:
|    o Pointer to an alnMatrix or alnMatrixTwoBit
|      structure with scores to filter.
| Output:
|  - Modifies:
|    o scoreAryL, startIndexAryUL, and endIndexAryUL to
|      have overlapping scores removed
\-------------------------------------------------------*/
#define altScoreFiltRef(mtrxSTPtr){\
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: altScoreFlitRef
   '  - Filters by removing alternative alignments that
   '    have an overlap in their position on reference.
   '    The overlap with the highest score is kept.
   '  o fun-07 sec-01:
   '    - Variable declerations
   '  o fun-07 sec-02:
   '    - Remove overlapping alignments
   '  o fun-07 sec-03:
   '    - Blank the empty/discared alignments
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   ulong numElmUL = (mtrxSTPtr)->lenArraysUL;\
   ulong lenRefUL = (mtrxSTPtr)->lenRefUL;\
   \
   long *scorePtrL = (mtrxSTPtr)->scoreAryL;\
   long *lastScorePL = scorePtrL + numElmUL - 1;\
   ulong *startPtrUL = (mtrxSTPtr)->startIndexAryUL;\
   ulong *endPtrUL = (mtrxSTPtr)->endIndexAryUL;\
   \
   long *scoreOnPL = (mtrxSTPtr)->scoreAryL;\
   ulong *startOnPUL = (mtrxSTPtr)->startIndexAryUL;\
   ulong *endOnPUL = (mtrxSTPtr)->endIndexAryUL;\
   \
   ulong nextRefStartUL = 0;\
   ulong refEndUL = 0;\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-02:
   ^  - Remove overlapping alignments
   ^  o fun-07 sec-02 sub-01:
   ^    - Loop setup + deal with blank or discared scores
   ^  o fun-07 sec-02 sub-02:
   ^    - Check if have an overlap
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*****************************************************\
   * Fun-07 Sec-02 Sub-01:
   *   - Loop setup + deal with blank or discared scores
   \****************************************************/\
   \
   sortRefScores((mtrxSTPtr), 0, numElmUL - 1);\
   \
   while(scoreOnPL <= lastScorePL)\
   { /*Loop: Check all scores to see if overlap*/\
      if(*scoreOnPL == 0)\
      { /*If: This had a score of 0 (was blank)*/\
         ++scoreOnPL;\
         ++startOnPUL;\
         ++endOnPUL;\
         continue;\
      } /*If: This had a score of 0 (was blank)*/\
      \
      if(*scorePtrL == 0)\
      { /*If: My comparision score is a blank*/\
         *scorePtrL = *scoreOnPL;\
         *startPtrUL = *startOnPUL;\
         *endPtrUL = *endOnPUL;\
         \
         ++scoreOnPL;\
         ++startOnPUL;\
         ++endOnPUL;\
         continue;\
      } /*If: My comparision score is a blank*/\
      \
      /**************************************************\
      * Fun-07 Sec-02 Sub-02:
      *   - Check if have an overlap
      \*************************************************/\
      \
      refEndUL = indexToRef(lenRefUL, *endPtrUL);\
      nextRefStartUL = indexToRef(lenRefUL, *startOnPUL);\
      \
      if(nextRefStartUL <= refEndUL)\
      { /*If: I have an overlap*/\
         if(*scorePtrL < *scoreOnPL)\
         { /*If: the new score is better*/\
            *scorePtrL = *scoreOnPL;\
            *startPtrUL = *startOnPUL;\
            *endPtrUL = *endOnPUL;\
         } /*If: the new score is better*/\
      } /*If: I have an overlap*/\
      \
      else if(scorePtrL != scoreOnPL)\
      { /*If: I am overwriting blank/overlapping score*/\
         ++scorePtrL;\
         ++startPtrUL;\
         ++endPtrUL;\
         *scorePtrL = *scoreOnPL;\
         *startPtrUL = *startOnPUL;\
         *endPtrUL = *endOnPUL;\
      } /*If: I am overwriting blank/overlapping score*/\
      \
      else\
      { /*Else: I need to move to the next score*/\
         ++scorePtrL;\
         ++startPtrUL;\
         ++endPtrUL;\
         /*case scorePtrL == scoreOnPL (same address)*/\
      } /*Else: I need to move to the next score*/\
      \
      ++scoreOnPL;\
      ++startOnPUL;\
      ++endOnPUL;\
   } /*Loop: Check all scores to see if overlap*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-03:
   ^  - Blank the empty/discared alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   if(scorePtrL < scoreOnPL)\
   { /*If: overwriting scores, get off last kept score*/\
      ++scorePtrL;\
      ++startPtrUL;\
      ++endPtrUL;\
   } /*If: overwriting scores, get off last kept score*/\
   \
   while(scorePtrL < scoreOnPL)\
   { /*Loop: blank ending scores*/\
      *scorePtrL = 0;\
      ++scorePtrL;\
   } /*Loop: blank ending scores*/\
} /*altScoreFiltRef*/

/*-------------------------------------------------------\
| Fun-08: altScoreFlitRefQry
|  - Filters by removing alternative alignments that
|    have an overlap in their position on the reference
|    or query. The overlap with the highest score is kept.
| Input:
|  - mtrxSTPtr:
|    o Pointer to an alnMatrix or alnMatrixTwoBit
|      structure with scores to filter.
| Output:
|  - Modifies:
|    o scoreAryL, startIndexAryUL, and endIndexAryUL to
|      have overlapping scores removed
\-------------------------------------------------------*/
#define altScoreFiltRefQry(mtrxSTPtr){\
   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: altScoreFlitRefQry
   '  - Filters by removing alternative alignments that
   '    have an overlap in their position on reference
   '    or the query. The overlap with the highest score
   '    is kept.
   '  o fun-08 sec-01:
   '    - Variable declerations
   '  o fun-08 sec-02:
   '    - Remove overlapping alignments
   '  o fun-08 sec-03:
   '    - Blank the empty/discared alignments
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-08 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   ulong numElmUL = (mtrxSTPtr)->lenArraysUL;\
   ulong lenRefUL = (mtrxSTPtr)->lenRefUL;\
   \
   long *scorePtrL = (mtrxSTPtr)->scoreAryL;\
   long *lastScorePL = scorePtrL + numElmUL - 1;\
   ulong *startPtrUL = (mtrxSTPtr)->startIndexAryUL;\
   ulong *endPtrUL = (mtrxSTPtr)->endIndexAryUL;\
   \
   long *scoreOnPL = (mtrxSTPtr)->scoreAryL;\
   ulong *startOnPUL = (mtrxSTPtr)->startIndexAryUL;\
   ulong *endOnPUL = (mtrxSTPtr)->endIndexAryUL;\
   \
   ulong nextRefStartUL = 0;\
   ulong refStartUL = 0;\
   ulong refEndUL = 0;\
   \
   ulong nextQryStartUL = 0;\
   ulong qryStartUL = 0;\
   ulong qryEndUL = 0;\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-08 Sec-02:
   ^  - Remove overlapping alignments
   ^  o fun-08 sec-02 sub-01:
   ^    - Loop setup + deal with blank or discared scores
   ^  o fun-08 sec-02 sub-02:
   ^    - Check if have an overlap
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   /*****************************************************\
   * Fun-08 Sec-02 Sub-01:
   *   - Loop setup + deal with blank or discared scores
   \****************************************************/\
   \
   sortRefQryScores((mtrxSTPtr), 0, numElmUL - 1);\
   \
   while(scoreOnPL <= lastScorePL)\
   { /*Loop: Check all scores to see if overlap*/\
      if(*scoreOnPL == 0)\
      { /*If: This had a score of 0 (was blank)*/\
         ++scoreOnPL;\
         ++startOnPUL;\
         ++endOnPUL;\
         continue;\
      } /*If: This had a score of 0 (was blank)*/\
      \
      if(*scorePtrL == 0)\
      { /*If: My comparision score is a blank*/\
         *scorePtrL = *scoreOnPL;\
         *startPtrUL = *startOnPUL;\
         *endPtrUL = *endOnPUL;\
         \
         ++scoreOnPL;\
         ++startOnPUL;\
         ++endOnPUL;\
         continue;\
      } /*If: My comparision score is a blank*/\
      \
      /**************************************************\
      * Fun-08 Sec-02 Sub-02:
      *   - Check if have an overlap
      \*************************************************/\
      \
      refStartUL = indexToRef(lenRefUL, *startPtrUL);\
      refEndUL = indexToRef(lenRefUL, *endPtrUL);\
      nextRefStartUL = indexToRef(lenRefUL, *startOnPUL);\
      \
      qryStartUL = indexToQry(lenRefUL, *startPtrUL);\
      qryEndUL = indexToQry(lenRefUL, *endPtrUL);\
      nextQryStartUL = indexToQry(lenRefUL, *startOnPUL);\
      \
      if(   nextRefStartUL <= refEndUL\
         && nextRefStartUL >= refStartUL\
         && nextQryStartUL <= qryEndUL\
         && nextQryStartUL >= qryStartUL\
      ){ /*If: I have an overlap*/\
         if(*scorePtrL < *scoreOnPL)\
         { /*If: the new score is better*/\
            *scorePtrL = *scoreOnPL;\
            *startPtrUL = *startOnPUL;\
            *endPtrUL = *endOnPUL;\
         } /*If: the new score is better*/\
      } /*If: I have an overlap*/\
      \
      else if(scorePtrL != scoreOnPL)\
      { /*If: I am overwriting blank/overlapping score*/\
         ++scorePtrL;\
         ++startPtrUL;\
         ++endPtrUL;\
         *scorePtrL = *scoreOnPL;\
         *startPtrUL = *startOnPUL;\
         *endPtrUL = *endOnPUL;\
      } /*If: I am overwriting blank/overlapping score*/\
      \
      else\
      { /*Else: I need to move to the next score*/\
         ++scorePtrL;\
         ++startPtrUL;\
         ++endPtrUL;\
         /*case scorePtrL == scoreOnPL (same address)*/\
      } /*Else: I need to move to the next score*/\
      \
      ++scoreOnPL;\
      ++startOnPUL;\
      ++endOnPUL;\
   } /*Loop: Check all scores to see if overlap*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-08 Sec-03:
   ^  - Blank the empty/discared alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   if(scorePtrL < scoreOnPL)\
   { /*If: overwriting scores, get off last kept score*/\
      ++scorePtrL;\
      ++startPtrUL;\
      ++endPtrUL;\
   } /*If: overwriting scores, get off last kept score*/\
   \
   while(scorePtrL < scoreOnPL)\
   { /*Loop: blank ending scores*/\
      *scorePtrL = 0;\
      ++scorePtrL;\
   } /*Loop: blank ending scores*/\
} /*altScoreFiltRefQry*/

#endif
