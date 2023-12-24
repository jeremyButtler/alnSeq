/*########################################################
# Name: twoBitArrays
# Use:
#   o Holds functions to handle two bit arrays
# Libraries:
#   - "dataTypeShortHand.h"
# C Standard Libraries:
#   - <stdlib.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - Has includes and defininitions
'  o st-01: twoBitArry
'    - Holds the two bit arrays and pointers to manage it
'  o fun-01: getCNegBit
'    - Get the negative flag for a character
'  o fun-02 getTwoBitElm:
'    - Get an element from a two bit array
'  o fun-02 twoBitAryShiftBytsForNewElm:
'    - Make room in unit8_t for two more two-bit elements
'  o fun-04 twoBitMvToNextElm:
'    - Moves to the next element in a two-bit array
'  o fun-05 twoBitMvForXElm:
'    - Moves two bit array pointer by x to next element
'  o fun-06 twoBitMvBackOneElm:
'    - Moves back one element in a 2-bit array
'  o fun-07 twoBitMvBackXElm:
'    - Moves back x elements in a 2-bit array
'  o fun-08 changeTwoBitElm:
'    - Changes a single two bit value in a two bit array.
'  o fun-09 blankTwoBitLimb:
'    - Sets a uint8_t (a limb) in a two bit array to 0.
'      Each limb holds four two bit elments.
'  o fun-10 twoBitMvToNextLimb:
'    - Moves to start of the next limb (uint8_t) in two
'      bit array
'  o fun-11 twoBitMvToLastLimb:
'    - Moves to start of the previous limb (uint8_t) in
'      two bit array
'  o fun-12 makeTwoBit:
'    - Make a two bit array struct
'  o fun-13 cpTwoBitPos:
'    - copy pointers of cpTwoBitST to dupTwoBitST
'  o fun-14 twoBitMvXElmFromStart:
'    - Change the two bit array postion to x elements
'       from the starting position
'  o fun-15 twoBitGetLen:
'    - Get the length of a two-bit array
'  o fun-16 twoBitGetIndex:
'    - Returns the index of the current two bit element
'  o fun-17 freeTwoBitStack:
'    - Frees the variables (array) in a twoBit structure
'  o fun-18 freeTwoBit:
'    - Frees a two bit struct and all variables in it
'  o fun-19 freeTwoBitStructOnly:
'    - Frees a two bit struct, but not the internal
'      variables, such as the two bit array
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Has includes and defininitions
\-------------------------------------------------------*/

#ifndef TWOBITARRAYS_H
#define TWOBITARRAYS_H

#include <stdlib.h>
#include "dataTypeShortHand.h"

/*-------------------------------------------------------\
| St-01: twoBitArry
|  - Holds the two bit arrays and pointers to manage it
\-------------------------------------------------------*/
typedef struct twoBitAry
{ /*twoBitArry structure*/
  uchar *firstLimbUCPtr; /*First limb in two bit array*/
  uchar *limbOnUCPtr;    /*Limb currently working on*/
  char elmOnC;           /*element on in the limb*/
  ulong lenAryUL;        /*Number of limbs in the array*/
}twoBitAry;

/*-------------------------------------------------------\
| Fun-01: getCNegBit
|  - Get the negative flag for a character
| Use:
|  - This is just to avoid (negChar < 0), which can have
|    branching issues on some cpus (likely very old
|    machines). It takes the same amount of time. The
|    (sizeof(unsigned char) << 3) - 1 statement will be
|    converted to a constant during compile time.
| Output:
|  - Returns:
|    o 1: if negChar was negative
|    o 0: if negChar was positive
\-------------------------------------------------------*/
#define getCNegBit(negChar)(\
   (( (uchar) (negChar) ) >> (((sizeof(uchar) << 3) - 1)))\
)

/*-------------------------------------------------------\
| Fun-02: getElmFromToBitUCAry
|  - Get an element from a two bit array
| Input:
|  - twoBitST
|     o Pointer to twoBitAry structure to get two bit
|       element from
| Output:
|  - Returns:
|    o Two bits of interest from the two bit array
| Note:
|  - Each limb has four two bit elements
\-------------------------------------------------------*/
#define getTwoBitElm(twoBitST)({\
    /*This was a branchless attempt. It had little affect
      on the speed. For this to work, the case must
      be reversed
    */\
    char shiftC = (twoBitST)->elmOnC << 1;\
    (*(twoBitST)->limbOnUCPtr & (3 << shiftC)) >> shiftC;\
}) /*getTwoBitElm*/

/*-------------------------------------------------------\
| Fun-04: twoBitMvToNextElm
|  - Moves to the next element in a two-bit array
| Input:
|  - twoBitST
|     o Pointer to twoBitAry structure to move to next
|       element in
| Output:
|  - Modifies:
|    o twoBitST to point to next element in two-bit array
\-------------------------------------------------------*/
#define twoBitMvToNextElm(twoBitST){\
   ++(twoBitST)->elmOnC;\
   (twoBitST)->limbOnUCPtr += ((twoBitST)->elmOnC >> 2);\
   (twoBitST)->elmOnC &= 3;\
} /*twoBitMvToNextElm*/

/*-------------------------------------------------------\
| Fun-05: twoBitMvForXElm
|  - Moves forward in two bit array by x elements
| Input:
|  - twoBitST
|     o Pointer to twoBitAry structure to move foward by
|       x elements
|  - shiftBy:
|     o How many element to move forward by (unsigned)
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements ahead
\-------------------------------------------------------*/
#define twoBitMvForXElm(twoBitST, shiftBy){\
   /*Get number of whole shifts to perform*/\
   (shiftBy) += (twoBitST)->elmOnC;\
   \
   /*Do the full shifts (by bytes)*/\
   (twoBitST)->limbOnUCPtr += ((shiftBy) >> 2);\
   \
   /*Find the element on in the byte (finish shift)*/\
   (twoBitST)->elmOnC = ((shiftBy) & 3);\
} /*twoBitMvForXElm*/

/*-------------------------------------------------------\
| Fun-06: twoBitMvBackOneElm
|  - Moves back one element in a 2-bit array
| Input:
|  - twoBitST
|     o Pointer to twoBitAry structure to move back one
|       element in
| Output:
|  - Modifies:
|    o twoBitST to point to the previous element
\-------------------------------------------------------*/
#define twoBitMvBackOneElm(twoBitST){\
   --(twoBitST)->elmOnC;\
   (twoBitST)->limbOnUCPtr -=\
      getCNegBit((twoBitST)->elmOnC);\
     /*Adds extra 1 if I am moving back extra elements*/\
   (twoBitST)->elmOnC = 3 & (twoBitST)->elmOnC;\
    /*-1 goes to 3, and all other values are 3 or less and
    ` so remain the same
    */\
} /*twoBitMvToNextElm*/

/*-------------------------------------------------------\
| Fun-07: twoBitMvBackXElm
|  - Moves back X elements back in a 2-bit array
| Input:
|  - twoBitST
|     o Pointer to twoBitAry structure to move backwards by
|       x elements
|  - shiftBy:
|     o How many element to move backwards by (unsigned)
| Output:
|  - Modifies:
|    o twoBotST to point to X elements back    
\-------------------------------------------------------*/
#define twoBitMvBackXElm(twoBitST, shiftBy){\
   /*Figure out how many elements I am moving back*/\
   (twoBitST)->elmOnC =\
      (twoBitST)->elmOnC - ((shiftBy) & 3);\
   \
   (twoBitST)->limbOnUCPtr -=\
      ((shiftBy) >> 2) + getCNegBit((twoBitST)->elmOnC);\
      /*Adds extra 1 if I am moving back extra elements*/\
   \
   (twoBitST)->elmOnC =\
        (twoBitST)->elmOnC\
      ^ (  -(getCNegBit((twoBitST)->elmOnC))\
         & ((twoBitST)->elmOnC ^ ((twoBitST)->elmOnC +4))\
        );\
     /*This is a minimize function, but instead of
     ` Returning the minimum if negative, it returns
     ` the minimum + 4.
     ` -(twoBitST->elmOnC < 0)
     `   Is -1 (111....) if elmOnC is negative
     `   Is 0 other if elmOnC is not negative.
     ` elmOnC + 4
     `   This is the correct ending index if elmOnC is
     `   negative.
     ` elmOnC ^ (-1 & (elmOnC ^ (elmOnC + 4)))
     `   This returns elmOnC + 4
     ` elmOnC ^ (0 & (elmOnC ^ (elmOnC + 4)))
     `   This returns elmOnC
     */\
} /*twoBitMvBackXElm*/

/*-------------------------------------------------------\
| Fun-08: changeTwoBitElm
|  - Changes a single two bit value in a two bit array.
| Input:
|  - twoBitST
|    o Pointer to twoBitAry structure to change
|  - element:
|    o element to change to (only first 2 bits used)
| Output:
|  - Modifies:
|    o Changes the taget element in twoBitUCArray to
|      input value
\-------------------------------------------------------*/
#define changeTwoBitElm(twoBitST, element){\
    char valC = (twoBitST)->elmOnC << 1;\
    char clearC = 3 << valC;\
    \
    valC = (element) << valC;\
    /*This shifts the new value to its correct position
    ` in the two bit array
    */\
    \
    *(twoBitST)->limbOnUCPtr =\
      (*(twoBitST)->limbOnUCPtr & (~clearC)) | valC;\
    /*This clears and then sets the value to its correct
    `  postion:
    ` & (~valC):
    `   ~valC removes
    ` and then
    ` sets the postion to the input value (| valC).
    */\
} /*changeTwoBitElm*/

/*-------------------------------------------------------\
| Fun-09: blankTwoBitLimb
|  - Sets a limb in a two bit array to 0.
|    Each limb holds four two bit elments.
| Input:
|  - twoBitST
|    o Pointer to twoBitAry structure with limb to blank
| Output:
|  - Modifies:
|    - Sets the current limb in a two-bit array to 0
\-------------------------------------------------------*/
#define blankTwoBitLimb(twoBitST){\
   (twoBitST)->limbOnUCPtr = 0;\
   (twoBitST)->elmOnC = 0;\
} /*blankTwoBitLimb*/

/*-------------------------------------------------------\
| Fun-10: twoBitMvToNextLimb
|  - Moves back to start of limb in a two bit array
| Input:
|  - twoBitST
|    o Pointer to twoBitAry structure to move to next limb
| Output:
|  - Modifies:
|    o Moves to the next limb in a two bit array
|  - Returns:
|    o 0: For a succesfull move
|    o 1: For memory out of bounds error
\-------------------------------------------------------*/
#define twoBitMvToNextLimb(twoBitST){\
   ++(twoBitST)->limbOnUCPtr;\
   (twoBitST)->elmOnC = 0;\
} /*twoBitMvToNextLimb*/

/*-------------------------------------------------------\
| Fun-11: twoBitMvToLastLimb
|  - Moves to the start of the previous limb (uint8_t)
|    in the two bit array
| Input:
|  - twoBitST
|    o Pointer to twoBitAry structure to move to last limb
| Output:
|  - Modifies:
|    o lastLimbUCPtr to point to the previous limb
\-------------------------------------------------------*/
#define twoBitMvToLastLimb(twoBitST){\
   --(twoBitST)->limbOnUCPtr;\
   (twoBitST)->elmOnC = 0;\
} /*twoBitMvToLastLimb*/

/*-------------------------------------------------------\
| Fun-12: makeTwoBit
|  - Make a two bit array struct
| Input:
|  - lenAry
|    o Length to make the two bit array
|  - blankAryBl
|    o 1: make a blank two bit structer
|    o 0: make an two bit strucuter with an array
| Output:
|  - Returns:
|    - twoBitAry structure with an unit8_t array of limbs
|    - blankAryBl = 1: returns a blank twoBitAry structer
\-------------------------------------------------------*/
#define makeTwoBit(lenAry, blankAryBl)({\
   struct twoBitAry *twoBitST = 0;\
   twoBitST = malloc(sizeof(struct twoBitAry));\
   \
   if(twoBitST != 0)\
   { /*If: I did not have an memroy error*/\
      twoBitST->elmOnC = 0;\
      \
      if(! (blankAryBl))\
      { /*If: I am making an array*/\
         twoBitST->firstLimbUCPtr =\
           calloc(((lenAry) >> 2) + 1, sizeof(uchar));\
           /*Each limb (uint8_t) has four elements, so I
           `  need to divide by 4 (>> 2)
           */\
         \
         if(twoBitST->firstLimbUCPtr == 0)\
         { /*If: I had a memory error*/\
            free(twoBitST);\
            twoBitST = 0;\
         } /*If: I had a memory error*/\
         \
         twoBitST->lenAryUL = (((lenAry) >> 2) << 2) + 4;\
         twoBitST->limbOnUCPtr =twoBitST->firstLimbUCPtr;\
      } /*If: I am making an array*/\
      \
      else\
      { /*Else: I am returing a blank structer*/\
         twoBitST->firstLimbUCPtr = 0;\
         twoBitST->limbOnUCPtr = 0;\
         twoBitST->lenAryUL = 0;\
      } /*Else: I am returing a blank structer*/\
   } /*If: I did not have an memroy error*/\
   \
   twoBitST; /*Is struct pointer or 0 for memory error*/\
}) /*makeTwoBitArray*/

/*-------------------------------------------------------\
| Fun-13: cpTwoBitPos
|  - copy pointers of cpTwoBitST to dupTwoBitST
| Input:
|  - cpTwoBitST:
|    o twoBit array to copy pointers from
|  - dupTwoBitST:
|    o twoBit array to copy pointerst to
| Output:
|  - Modifies:
|    o dupTwoBitST to hold same pointers/values as
|      cpTwoBitST
\-------------------------------------------------------*/
#define cpTwoBitPos(cpTwoBitST, dupTwoBitST){\
   (dupTwoBitST)->firstLimbUCPtr =\
      (cpTwoBitST)->firstLimbUCPtr;\
   \
   (dupTwoBitST)->limbOnUCPtr = (cpTwoBitST)->limbOnUCPtr;\
   (dupTwoBitST)->elmOnC = (cpTwoBitST)->elmOnC;\
   (dupTwoBitST)->lenAryUL = (cpTwoBitST)->lenAryUL;\
} /*cpTwoBitPos*/

/*-------------------------------------------------------\
| Fun-14: twoBitMvXElmFromStart
|  - Change the two bit array postion to x elements
|    from the starting position
| Input:
|  - twoBitST:
|    o pointer to twoBit array to change position in
|  - shiftBy:
|    o Number of elements to jump to from start
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements after the start of
|      the array
\-------------------------------------------------------*/
#define twoBitMvXElmFromStart(twoBitST, shiftBy){\
   /*Get the number of whole shifts to do*/\
   (twoBitST)->limbOnUCPtr =\
     (twoBitST)->firstLimbUCPtr + ((shiftBy) >> 2);\
   \
   (twoBitST)->elmOnC = ((shiftBy) & 3);\
} /*twoBitMvXElmFromStart*/

/*-------------------------------------------------------\
| Fun-15: twoBitGetLen
|  - Get the length of a two-bit array
| Input:
|  - twoBitST:
|    o pointer to twoBit array to get length for
| Output:
|  - Returns:
|    o Length of the two bit array
\-------------------------------------------------------*/
#define twoBitGetLen(twoBitST) ((twoBitST)->lenAryUL << 2)

/*-------------------------------------------------------\
| Fun-16: twoBitGetIndex
|  - Returns the index of the current two bit element
|    on
| Input:
|  - twoBitST:
|    o pointer to twoBit array to get length for
| Output:
|  - Returns
|    o The index of the elemen on in the twoBitST struct
\-------------------------------------------------------*/
#define twoBitGetIndex(twoBitST)(\
  (\
    ( (\
          (twoBitST)->limbOnUCPtr\
        - (twoBitST)->firstLimbUCPtr\
      ) << 2\
    )\
   + (twoBitST)->elmOnC\
  ) /*Return the index*/\
) /*twoBitGetIndex*/

/*-------------------------------------------------------\
| Fun-17: freeTwoBitStack
|  - Frees the interal variables in a two bit array
| Input:
|  - twoBitSTPtr:
|    o pointer to twoBit array with variables to free
| Output:
|  - Frees:
|    o array in twoBitST
|  - Sets:
|    o All values in twoBitST to 0
\-------------------------------------------------------*/
#define freeTwoBitStack(twoBitSTPtr){\
   if((twoBitSTPtr) != 0)\
      if((twoBitSTPtr)->firstLimbUCPtr != 0)\
      { /*If: I am freeing the array*/\
         free((twoBitSTPtr)->firstLimbUCPtr);\
         (twoBitSTPtr)->firstLimbUCPtr = 0;\
         (twoBitSTPtr)->limbOnUCPtr = 0;\
         (twoBitSTPtr)->elmOnC = 0;\
         (twoBitSTPtr)->lenAryUL = 0;\
      } /*If: I am freeing the array*/\
} /*freeTwoBitStack*/

/*-------------------------------------------------------\
| Fun-18: freeTwoBit
|  - Frees a two bit array
| Input:
|  - twoBitSTPtr:
|    o pointer to twoBit structure to free
| Output:
|  - Frees:
|    o twoBitSTPtr and all variables in twoBitSTPtr
|  - Sets:
|    o twoBitSTPtr to 0
\-------------------------------------------------------*/
#define freeTwoBit(twoBitSTPtr){\
   if((twoBitSTPtr) != 0)\
   { /*If: I need to free twoBitSTPtr*/\
      freeTwoBitStack((twoBitSTPtr))\
      free(twoBitSTPtr);\
      (twoBitSTPtr) = 0;\
   } /*If: I need to free twoBitSTPtr*/\
} /*freeTwoBit*/

/*-------------------------------------------------------\
| Fun-19: freeTwoBitStructOnly
|  - Frees the two bit structer, but no internal variables
| Input:
|  - twoBitSTPtr:
|    o pointer to twoBit structure to free
| Output:
|  - Frees:
|    o twoBitSTPtr structure, but not the twoBit array
|  - Sets:
|    o twoBitSTPtr to 0
\-------------------------------------------------------*/
#define freeTwoBitStructOnly(twoBitSTPtr){\
   if((twoBitSTPtr) != 0)\
   { /*If: I need to free twoBitSTPtr*/\
      free(twoBitSTPtr);\
      (twoBitSTPtr) = 0;\
   } /*If: I need to free twoBitSTPtr*/\
} /*freeTwoBitStructOnly*/

#endif
