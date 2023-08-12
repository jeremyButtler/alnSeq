/*#########################################################
# Name: twoBitArrays
# Use:
#   o Holds functions to handle two bit arrays
# Includes:
# C Standard Includes:
#   - <stdint.h>
#   - <stdlib.h>
#########################################################*/

#ifndef TWOBITARRAYS_H
#define TWOBITARRAYS_H

#include <stdint.h>
#include <stdlib.h>

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  - struct-01 twoBitArry:
'     o Holds the two bit arrays and pointers to manage it
'  - fun-01 getElmFromToBitUCAry:
'     o Get an element from a two bit array
'  - fun-03 twoBitAryMoveToNextElm:
'     o Moves to the next element in a two-bit array
'  - fun-04 twoBitAryMoveForXElm:
'     o Moves two bit array pointer by x to next element
'  - fun-05 twoBitAryMoveBackOneElm:
'     o Moves back one element in a 2-bit array
'  - fun-06 twoBitAryMoveBackXElm:
'     o Moves back x elements in a 2-bit array
'  - fun-08 changeTwoBitElm:
'     o Changes a single two bit value in a two bit array.
'  - fun-09 blankLimb:
'     o Sets a uint8_t (a limb) in a two bit array to 0.
'       Each limb holds four two bit elments.
'  - fun-10 moveToNextLimb:
'     o Moves to start of the next limb (uint8_t) in two
'       bit array
'  - fun-11 moveToLastLimb:
'     o Moves to start of the previous limb (uint8_t) in
'       two bit array
'  - fun-12 makeTwoBitArry:
'     o Make a two bit array struct
'  - fun-13 cpTwoBitPos:
'     o copy pointers of cpTwoBitST to dupTwoBitST
'  - fun-14 freeTwoBitAry:
'     o Frees a two bit array
'  - fun-15 moveXElmFromStart:
'     o Change the two bit array postion to x elements
'       from the starting position
'  - fun-16 lenTwoBitAry:
'     o Get the length of a two-bit array
'  - Fun-17 TOC: Sec-01: getTwoBitAryIndex
'    - Returns the index of the current two bit element
'
' Note: The remaing functions are functions that move
'       foward and backwards and check if the movement is
'       out of bounds. The are copies of previous functions
'
'  - fun-18 twoBitAryMoveToNextElmBoundsCheck:
'     o Moves to the next element in a two-bit array and
'       returns error if out of bounds
'  - fun-19 twoBitAryMoveForXElmBoundsCheck:
'     o Moves forward in two bit array by x elements
'  - fun-20 twoBitAryMoveBackOneElm
'     o Moves back one element in a 2-bit array
'  - fun-21 twoBitAryMoveBackXElmBoundsCheck:
'    - Moves back X elements back in a 2-bit array
'  - fun-22 Sub-1: moveToNextLimbBoundsCheck:
'     o Moves back to start of limb in a two bit array
'  - fun-23 moveToLastLimbBoundsCheck:
'     o Moves to the start of the previous limb (uint8_t)
'       in the two bit array
'  - fun-24 moveXElmFromStartBoundsCheck:
'     o Change the two bit array postion to x elements
'       from the starting position
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Interfacing: A quick guid to setting up a two bit array
|   - This is here to describe the order you might want to
|     use these functions.
|   o Make an uint8_t array for the twobit array
|     - Size should be: (maximum elemnents / 4) + 1
|     - This is an array of limbs, were each limb holds
|       four elements
|   o make a uint8_t bit counter set to 0
|     (uint8_t bitUC = 0)
|     - This is to keep track of the element on in the limb
|   o TwoBitAryShiftBytsForNewElm is here to allow you to
|     manipulate two bit array elements directly. You can
|     set new elements by *twoBitArray |= newTwoBitElement
|     and then calling TwoBitAryShiftBytsForNewElm to
|     shift the elements.
|     - This will result in the very end elements of the
|       array not being shifted correctly, so you will want
|       to call finishTowBitAry at the end.
|   o An alternative way to add elements would be call
|     changeTwoBitElm, which replaces the old elements with
|     the target elements.
|   o After adding an element you will want to move to the
|     next element with twoBitAryMoveToNextElm
|     - Use twoBitAryMoveBackOneElm to move backwards
|   o Your can check an element using getElmFromToBitUCAry
\--------------------------------------------------------*/

/*--------------------------------------------------------\
| Struct-01: twoBitArry
|  - Holds the two bit arrays and pointers to manage it
\--------------------------------------------------------*/
typedef struct twoBitAry
{ // twoBitArry structure
  uint8_t *firstLimbUCPtr; // First limb in two bit array
  uint8_t *limbOnUCPtr;    // Limb currently working on
  uint8_t elmOnUC;         // element on in the limb
  unsigned long lenAryUL;  // Number of limbs in the array
}twoBitAry;

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Two bits of interest from the two bit array
| Note:
|  - Each limb has four two bit elements
\--------------------------------------------------------*/
uint8_t getTwoBitAryElm(
  struct twoBitAry *twoBitST // Array to get element from
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: getElmFromToBitUCAry
   '  - Get an element from a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to next element in two-bit array
\--------------------------------------------------------*/
void twoBitAryMoveToNextElm(
  struct twoBitAry *twoBitST
    // Two bit array to change index
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: twoBitAryMoveToNextElm
   '  - Moves to the next element in a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements ahead
\--------------------------------------------------------*/
void twoBitAryMoveForXElm(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL    // Number elements to shift by
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: twoBitAryMoveForXElm
   '  - Moves forward in two bit array by x elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
     
/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to the previous element
\--------------------------------------------------------*/
void twoBitAryMoveBackOneElm(
  struct twoBitAry *twoBitST // array to move back in
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-1 Sub-1: twoBitAryMoveBackOneElm
   '  - Moves back one element in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBotST to point to X elements back    
\--------------------------------------------------------*/
void twoBitAryMoveBackXElm(
  struct twoBitAry *twoBitST, // To bit array to move back
  unsigned long shiftByUL // number elements to shift back
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: twoBitAryMoveBackXElm
   '  - Moves back X elements back in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Changes the taget element in twoBitUCArray to
|      input value
\--------------------------------------------------------*/
void changeTwoBitElm(
  struct twoBitAry *twoBitST, //Points to element to change
  uint8_t newValueUC          // New value; only the first
                              // two bits can be set to 1
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-1 Sub-1: changeTwoBitElm
   '  - Changes a single two bit value in a two bit array.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    - Sets the current limb in a two-bit array to 0
\--------------------------------------------------------*/
void blankLimb(
  struct twoBitAry *twoBitSt // two bit array to restart
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-1 Sub-1: blankLimb
   '  - Sets a limb in a two bit array to 0.
   '    Each limb holds four two bit elments.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Moves to the next limb in a two bit array
\--------------------------------------------------------*/
void moveToNextLimb(
  struct twoBitAry *twoBitST // two-bit array to move back
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-1 Sub-1: moveToNextLimb
   '  - Moves back to start of limb in a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o lastLimbUCPtr to point to the previous limb
\--------------------------------------------------------*/
void moveToLastLimb(
    struct twoBitAry *twoBitST // array to move back in
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: Sec-1 Sub-1: moveToLastLimb
   '  - Moves to the start of the previous limb (uint8_t)
   '     in the two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    - twoBitAry structure with an unit8_t array of limbs
|    - blankAryBl = 1: returns a blank twoBitAry structer
\--------------------------------------------------------*/
struct twoBitAry * makeTwoBitArry(
  unsigned long lenArryUL, // Length of new array
  char blankAryBl          // 1: Make a blank strucuter
                           // 0: Make structure with array
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: makeTwoBitArry
   '  - Make a two bit array struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o dupTwoBitST to hold same pointers/values as
|      cpTwoBitST
\--------------------------------------------------------*/
void cpTwoBitPos(
  struct twoBitAry *cpTwoBitST,  // Structer to copy
  struct twoBitAry *dupTwoBitST  // Struct to copy to
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: Sec-1 Sub-1: cpTwoBitPos
   '  - copy pointers of cpTwoBitST to dupTwoBitST
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  Frees: the input towBitAry strucuter
\--------------------------------------------------------*/
void freeTwoBitAry(
  struct twoBitAry *stToFree, // Two bit array to free
  char twoBitOnStackBl,       // 0: free stToFree
                              // 1: stToFree on stack
  char doNotFreeAryBl         // 0: Free everthing
      // 1: Do not free the structer in the array. Only
      //    use this if this stuct is a copy
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: freeTwoBitAry
   '  - Frees a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements after the start of
|      the array
\--------------------------------------------------------*/
void moveXElmFromStart(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL    // Number elements to shift by
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-15 TOC: Sec-1 Sub-1: moveXElmFromStart
   '  - Change the two bit array postion to x elements
   '    from the starting position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Length of the two bit array
\--------------------------------------------------------*/
unsigned long lenTwoBitAry(
  struct twoBitAry *twoBitST
    // two-bit array to get length for
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: lenTwoBitAry
   '  - Get the length of a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o The index of the elemen on in the twoBitST struct
\--------------------------------------------------------*/
unsigned long getTwoBitAryIndex(
  struct twoBitAry *twoBitST
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-17 TOC: Sec-01: getTwoBitAryIndex
   '  - Returns the index of the current two bit element
   '    on
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to next element in two-bit array
|  - Returns:
|    o 0: if could move to the next element
|    o 1: if could not move to next element
|      - on the last element in the array
\--------------------------------------------------------*/
char twoBitAryMoveToNextElmBoundsCheck(
  struct twoBitAry *twoBitST
    // Two bit array to change index
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-18 TOC: Sec-1: twoBitAryMoveToNextElmBoundsCheck
   '  - Moves to the next element in a two-bit array and
   '    returns error if out of bounds
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements ahead
|  - Returns:
|    o 0: For a sucessfull shift
|    o 1: For array out of bounds
\--------------------------------------------------------*/
char twoBitAryMoveForXElmBoundsCheck(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL    // Number elements to shift by
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-19 TOC: Sec-1 twoBitAryMoveForXElmBoundsCheck
   '  - Moves forward in two bit array by x elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to the previous element
|  - Returns:
|    o 0: For success
|    o 1: For array out of bounds
|    o 2: error; limb element is not 0-3 (elmOnUC)
\--------------------------------------------------------*/
char twoBitAryMoveBackOneElmBoundsCheck(
  struct twoBitAry *twoBitST // array to move back in
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-20 TOC: Sec-1 Sub-1: twoBitAryMoveBackOneElm
   '  - Moves back one element in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBotST to point to X elements back    
|  - Returns:
|    o 0: For success
|    o 1: For array out of bounds
|    o 2: error; limb element is not 0-3 (elmOnUC)
\--------------------------------------------------------*/
char twoBitAryMoveBackXElmBoundsCheck(
  struct twoBitAry *twoBitST, // To bit array to move back
  unsigned long shiftByUL  // number elements to shift back
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-21 TOC: Sec-1 twoBitAryMoveBackXElmBoundsCheck
   '  - Moves back X elements back in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Moves to the next limb in a two bit array
|  - Returns:
|    o 0: For a succesfull move
|    o 1: For memory out of bounds error
\--------------------------------------------------------*/
char moveToNextLimbBoundsCheck(
  struct twoBitAry *twoBitST // two-bit array to move back
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-22 TOC: Sec-1 Sub-1: moveToNextLimbBoundsCheck
   '  - Moves back to start of limb in a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o lastLimbUCPtr to point to the previous limb
|  - Returns:
|    o 0: if could move back on limb
|    o 1: if I am on the first element 
|    o 2: if I was on the frist limb, but not first element
\--------------------------------------------------------*/
char moveToLastLimbBoundsCheck(
    struct twoBitAry *twoBitST // array to move back in
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-23 TOC: Sec-1: moveToLastLimbBoundsCheck
   '  - Moves to the start of the previous limb (uint8_t)
   '    in the two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements after the start of
|      the array
|  - Returns:
|    o 0: For a sucessfull shift
|    o 1: For array out of bounds
\--------------------------------------------------------*/
char moveXElmFromStartBoundsCheck(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL  // Number elements to shift by
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-24 TOC: Sec-1: moveXElmFromStartBoundsCheck
   '  - Change the two bit array postion to x elements
   '    from the starting position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif

