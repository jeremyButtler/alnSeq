/*########################################################
# Name: alnMatrixStruct
#  - Holds the alnMatrix structures and functions
# Libraries:
#  - "twoBitArrays.h"
#  - "dataTypeShortHand.h"
# C Standard Libraries:
#  o <stdlib.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o st-01 alnMatrixStruct:
'    - Holds the direction matrix and best score(s) for a
'      single aligment
'  o st-02 alnMatixTwoBit:
'    - Holds the direction matrix and best score(s) for a
'      single aligment using a two bit matrix
'  o fun-01 initAlnMatrix:
'    - Sets all variables in matrixST to 0
'  o fun-02 initAlnMatrixTwoBit:
'    - Initializes an alnMatrixTwoBit struture
'  o fun-03 freeAlnMatrixStack:
'    - Frees all variables in an alnMatrix structure
'  o fun-04 freeAlnMatrixTwoBitStack:
'    - Frees all variables in an alnMatrixTwoBit structure
'  o fun-05 freeAlnMatrix:
'    - Frees an alnMatrix structure & its variables
'  o fun-06 freeAlnMatrixTwoBit:
'    - Frees an alnMatrixTwoBit structure & its variables
'  o fun-07 indexToQry:
'    - Gets the query coordinates of the query sequence in
'      an matrix.
'  o fun-08 indexToRef:
'    - Gets the coordinates of the reference sequence in
'      an matrix.
'  o fun-09 indexToCoord:
'    - Gets the coordinates of the reference and query
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef ALNMATRIXSTRUCT_H
#define ALNMATRIXSTRUCT_H

#include "twoBitArrays.h"
#include "dataTypeShortHand.h"

/*-------------------------------------------------------\
| St-01: alnMatrix
|  - Holds the direction matrix and best score(s) for a
|    single aligment
\-------------------------------------------------------*/
typedef struct alnMatrix
{ /*alnStruct*/
  char *dirMatrix;

  long bestScoreL;
  ulong bestStartIndexUL;
  ulong bestEndIndexUL;

  ulong lenRefUL;
  ulong refOffsetUL;

  ulong lenQryUL;
  ulong qryOffsetUL;

  ulong lenArraysUL;
  long *scoreAryL;
  ulong *startIndexAryUL;
  ulong *endIndexAryUL;
}alnMatrixStruct;

/*-------------------------------------------------------\
| St-02: alnMatixTwoBit
|  - Holds the direction matrix and best score(s) for a
|    single aligment using a two bit matrix
\-------------------------------------------------------*/
typedef struct alnMatrixTwoBit
{ /*alnStruct*/
  struct twoBitAry *dirMatrix;

  long bestScoreL;
  ulong bestStartIndexUL;
  ulong bestEndIndexUL;

  ulong lenRefUL;
  ulong refOffsetUL;

  ulong lenQryUL;
  ulong qryOffsetUL;

  ulong lenArraysUL;

  long *scoreAryL;
  ulong *startIndexAryUL;
  ulong *endIndexAryUL;
}alnMatrixStructTwoBit;

/*-------------------------------------------------------\
| Fun-01: initAlnMatrix
|  - Initializes an alnMatrix struture
| Input:
|  - matrixSTPtr:
|    o A pointer to an alnMatrix or alnMatrixTwoBit
|      structure initialize to
| Output:
|  - Sets:
|    o All variables in matrixSTPtr to 0
\-------------------------------------------------------*/
#define initAlnMatrix(matrixSTPtr){\
   (matrixSTPtr)->dirMatrix = 0;\
   \
   (matrixSTPtr)->bestScoreL = 0;\
   (matrixSTPtr)->bestEndIndexUL = 0;\
   (matrixSTPtr)->bestStartIndexUL = 0;\
   \
   (matrixSTPtr)->lenRefUL = 0;\
   (matrixSTPtr)->refOffsetUL = 0;\
   (matrixSTPtr)->lenQryUL = 0;\
   (matrixSTPtr)->qryOffsetUL = 0;\
   \
   (matrixSTPtr)->lenArraysUL = 0;\
   (matrixSTPtr)->scoreAryL = 0;\
   (matrixSTPtr)->startIndexAryUL = 0;\
   (matrixSTPtr)->endIndexAryUL = 0;\
} /*initAlnMatrixST*/

/*-------------------------------------------------------\
| Fun-02: initAlnMatrixTwoBit
|  - Initializes an alnMatrixTwoBit struture
| Input:
|  - matrixSTPtr:
|    o A pointer to an alnMatrixTwoBit structure to
|      initialize
| Output:
|  - Sets:
|    o All variables in matrixSTPtr to 0
\-------------------------------------------------------*/
#define initAlnMatrixTwoBit(matrixSTPtr){\
   (matrixSTPtr)->dirMatrix = 0;\
   \
   (matrixSTPtr)->bestScoreL = 0;\
   (matrixSTPtr)->bestEndIndexUL = 0;\
   (matrixSTPtr)->bestStartIndexUL = 0;\
   \
   (matrixSTPtr)->lenRefUL = 0;\
   (matrixSTPtr)->refOffsetUL = 0;\
   (matrixSTPtr)->lenQryUL = 0;\
   (matrixSTPtr)->qryOffsetUL = 0;\
   \
   (matrixSTPtr)->lenArraysUL = 0;\
   (matrixSTPtr)->scoreAryL = 0;\
   (matrixSTPtr)->startIndexAryUL = 0;\
   (matrixSTPtr)->endIndexAryUL = 0;\
} /*initAlnMatrixST*/

/*-------------------------------------------------------\
| Fun-03: freeAlnMatrixStack
|  - Frees all variable in an alnMatrix structure
| Input:
|  - matrixSTPtr
|    o Pointer to alnMatrix structure with variables to
|      free
| Output:
|  - Frees:
|    o All array/matrix variables in alnMatrixSTPtr
|  - Sets:
|    o All array/matrix variables in alnMatrixSTPtr to
|      0
\-------------------------------------------------------*/
#define freeAlnMatrixStack(matrixSTPtr){\
   if((matrixSTPtr) != 0)\
   { /*If: this is a structure (not null)*/\
      if((matrixSTPtr)->dirMatrix != 0)\
      { /*If: I need to free the direction matrix*/\
         free((matrixSTPtr)->dirMatrix);\
         (matrixSTPtr)->dirMatrix = 0;\
      } /*If: I need to free the direction matrix*/\
      \
      if((matrixSTPtr)->scoreAryL != 0)\
      { /*If: I need to free the reference scores*/\
         free((matrixSTPtr)->scoreAryL);\
         (matrixSTPtr)->scoreAryL = 0;\
      } /*If: I need to free the reference scores*/\
      \
      if((matrixSTPtr)->startIndexAryUL != 0)\
      { /*If: I need to free the reference start indexs*/\
         free((matrixSTPtr)->startIndexAryUL);\
         (matrixSTPtr)->startIndexAryUL = 0;\
      } /*If: I need to free the reference start indexs*/\
      \
      if((matrixSTPtr)->endIndexAryUL != 0)\
      { /*If: I need to free the qryerence end indexs*/\
         free((matrixSTPtr)->endIndexAryUL);\
         (matrixSTPtr)->endIndexAryUL = 0;\
      } /*If: I need to free the qryerence end indexs*/\
      \
      (matrixSTPtr)->bestScoreL = 0;\
      (matrixSTPtr)->bestEndIndexUL = 0;\
      (matrixSTPtr)->bestStartIndexUL = 0;\
      \
      (matrixSTPtr)->lenRefUL = 0;\
      (matrixSTPtr)->refOffsetUL = 0;\
      \
      (matrixSTPtr)->lenQryUL = 0;\
      (matrixSTPtr)->qryOffsetUL = 0;\
      \
      (matrixSTPtr)->lenArraysUL = 0;\
   } /*If: this is a structure (not null)*/\
} /*freeAlnMatrixStack*/

/*-------------------------------------------------------\
| Fun-04: freeAlnMatrixTwoBitStack
|  - Frees all variable in an alnMatrixTwo structure
| Input:
|  - matrixSTPtr
|    o Pointer to alnMatrix structure with variables to
|      free
| Output:
|  - Frees:
|    o All array/matrix variables in alnMatrixSTPtr
|  - Sets:
|    o All array/matrix variables in alnMatrixSTPtr to
|      0
\-------------------------------------------------------*/
#define freeAlnMatrixTwoBitStack(matrixSTPtr){\
   if((matrixSTPtr) != 0)\
   { /*If: this is a structure (not null)*/\
      if((matrixSTPtr)->dirMatrix != 0)\
         freeTwoBit((matrixSTPtr)->dirMatrix);\
      \
      if((matrixSTPtr)->scoreAryL != 0)\
      { /*If: I need to free the reference scores*/\
         free((matrixSTPtr)->scoreAryL);\
         (matrixSTPtr)->scoreAryL = 0;\
      } /*If: I need to free the reference scores*/\
      \
      if((matrixSTPtr)->startIndexAryUL != 0)\
      { /*If: I need to free the reference start indexs*/\
         free((matrixSTPtr)->startIndexAryUL);\
         (matrixSTPtr)->startIndexAryUL = 0;\
      } /*If: I need to free the reference start indexs*/\
      \
      if((matrixSTPtr)->endIndexAryUL != 0)\
      { /*If: I need to free the qryerence end indexs*/\
         free((matrixSTPtr)->endIndexAryUL);\
         (matrixSTPtr)->endIndexAryUL = 0;\
      } /*If: I need to free the qryerence end indexs*/\
      \
      (matrixSTPtr)->bestScoreL = 0;\
      (matrixSTPtr)->bestEndIndexUL = 0;\
      (matrixSTPtr)->bestStartIndexUL = 0;\
      \
      (matrixSTPtr)->lenRefUL = 0;\
      (matrixSTPtr)->refOffsetUL = 0;\
      \
      (matrixSTPtr)->lenQryUL = 0;\
      (matrixSTPtr)->qryOffsetUL = 0;\
      \
      (matrixSTPtr)->lenArraysUL = 0;\
   } /*If: this is a structure (not null)*/\
} /*freeAlnMatrixStack*/

/*-------------------------------------------------------/
| Fun-05: freeAlnMatrix
|  - Frees an alnMatrix structure and its variables
| Input:
|  - matrixSTPtr
|    o Pointer to alnMatrix structure to free
| Output:
|  - Frees:
|    o alnMatrixPtr
|  - Sets:
|    o alnMatrixPtr to 0
\-------------------------------------------------------*/
#define freeAlnMatrix(alnMatrixSTPtr){\
   freeAlnMatrixStack((alnMatrixSTPtr));\
   if((alnMatrixSTPtr) != 0) free((alnMatrixSTPtr));\
   (alnMatrixSTPtr) = 0;\
} /*freeAlnMatrix*/

/*-------------------------------------------------------/
| Fun-06: freeAlnMatrixTwoBit
|  - Frees an alnMatrixTwoBit structure and its variables
| Input:
|  - matrixSTPtr
|    o Pointer to alnMatrixTwoBit structure to free
| Output:
|  - Frees:
|    o alnMatrixPtr
|  - Sets:
|    o alnMatrixPtr to 0
\-------------------------------------------------------*/
#define freeAlnMatrixTwoBit(alnMatrixSTPtr){\
   freeAlnMatrixTwoBitStack((alnMatrixSTPtr));\
   if((alnMatrixSTPtr) != 0) free((alnMatrixSTPtr));\
   (alnMatrixSTPtr) = 0;\
} /*freeAlnMatrixTwoBit*/

/*-------------------------------------------------------\
| Fun-07: indexToQry
|   - Gets the query coordinates of the query sequence in
|     an matrix.
| Input
|   - refLen:
|     o Length of the reference sequence
|   - index:
|     o Index to convert to cooridnates
|   - qryCoord:
|     o Will hold the coordinate of the query sequence
| Output:
|  - Returns:
|    o The query coordinate that was in the index
\-------------------------------------------------------*/
#define indexToQry(refLen, index)(\
   ((index) / ((refLen) + 1)) - ((index) >= (refLen))\
   /* Find query coordinates from an index:
   `  - refLen + 1:
   `    o gives the length of each row in the matrix
   `    o The + 1 accounts for the gap column
   `  - index / (refLen + 1)
   `    o The number of rows down, which is the number of
   `      query bases + the gap row
   `  - index >= refLen
   `    o Is 1 if the index is divisible by the reference
   `      This means I have an index 1 value
   `    o else index is in the gap row (index 0)
   `  - position - (index 1 or 0)
   `    o Removes the gap row from the query if it is
   `      present. This makes sure I get an index 0 result
   */\
) /*indexToQry*/

/*-------------------------------------------------------\
| Fun-08: indexToRef
|   - Gets the coordinates of the reference sequence in
|     an matrix.
| Input
|   - refLen:
|     o Length of the reference sequence
|   - index:
|     o Index to convert to cooridnates
|   - refCoord:
|     o Will hold the coordinate of the reference sequence
| Output:
|  - Returns
|    o The reference coordinate that was in the index
\-------------------------------------------------------*/
#define indexToRef(refLen, index)(\
     ((index) % ((refLen) + 1))\
   - (((index) / ((refLen) + 1)) > 0)\
   /* Find reference coordinates from an index:
   `  - refLen + 1:
   `    o gives the length of each row in the matrix
   `    o The + 1 accounts for the gap column
   `  - index % (refLen + 1)
   `    o The number of columns, which is the number of
   `      reference bases + the gap row
   `  - index / (refLen + 1)
   `    o Gets the query position
   `  - (index / (refLen + 1)) > 0
   `    o Is 0 if I am in the gap column, 1 if on a real
   `      reference base
   `  - position - (index 1 or 0)
   `    o Removes the gap column from the reference if it
   `      is present. This makes sure I get an index 0
   `      result.
   */\
) /*indexToRef*/

/*-------------------------------------------------------\
| Fun-09: indexToCoord
|   - Gets the coordinates of the reference and query
|     sequence in an matrix.
| Input
|   - refLen:
|     o Length of the reference sequence
|   - index:
|     o Index to convert to cooridnates
|   - refCoord:
|     o Will hold the coordinate of the reference sequence
|   - qryCoord:
|     o Will hold the coordinate of the query sequence
| Output:
|  - Sets
|    o refCoord to the reference coordinate in index
|    o qryCoord to the query coordinate in index
\-------------------------------------------------------*/
#define indexToCoord(refLen, index, refCoord, qryCoord){\
   (refCoord) = indexToRef((refLen), (index));\
   (qryCoord) = indexToQry((refLen), (index));\
} /*indexToCoord*/


#endif
