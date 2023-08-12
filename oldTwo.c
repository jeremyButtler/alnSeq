/*#########################################################
# Name: twoBitArrays
# Use:
#   o Holds functions to handle two bit arrays
# Includes:
# C Standard Includes:
#   - <stdint.h>
#########################################################*/

#include "twoBitArrays.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
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
'  - fun-17 getTwoBitAryIndex:
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
| Output:
|  - Returns:
|    o Two bits of interest from the two bit array
| Note:
|  - Each limb has four two bit elements
\--------------------------------------------------------*/
uint8_t getTwoBitAryElm(
  struct twoBitAry *twoBitST /*Array to get element from*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: getElmFromToBitUCAry
   '  - Get an element from a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*This has a very small speed boost*/
    char shiftC = twoBitST->elmOnUC << 1;
    return (*twoBitST->limbOnUCPtr &(3<< shiftC)) >>shiftC;

    /* This Returns a chacater with all bits except the
    `    target bits cleared
    ` twoBitST->elmOnUC << 1:
    `   Gets position to shift bits to
    ` return statement:
    `   & 3 << shiftC:
    `     Clears bits not in target secton
    `   >> shiftC:
    `     Moves target bits to the zero position
    */
} /*getTwoBitAryElm*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to next element in two-bit array
\--------------------------------------------------------*/
void twoBitAryMoveToNextElm(
  struct twoBitAry *twoBitST
    // Two bit array to change index
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: twoBitAryMoveToNextElm
   '  - Moves to the next element in a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(twoBitST->elmOnUC)
   { // Switch; check if need to move to a new element
       case 0:
       case 1:
       case 2:
          ++twoBitST->elmOnUC; // Move to next element
          return;
       case 3:            // Next element is in next limb
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 0;
          return;
   } // Switch; check if need to move to a new element

   return;
} /*twoBitAryMoveToNextElm*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements ahead
\--------------------------------------------------------*/
void twoBitAryMoveForXElm(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL    // Number elements to shift by
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: twoBitAryMoveForXElm
   '  - Moves forward in two bit array by x elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Get whole shifts to perform
   twoBitST->limbOnUCPtr += (shiftByUL >> 2);
   //shiftByUL += twoBitST->elmOnUC;
   //twoBitST->limbOnUCPtr += (shiftByUL >> 2);
   //twoBitST->elmOnUC = shiftByUL & (1 | 2);
   //return;

   // Move the bit above 
   switch(twoBitST->elmOnUC + (shiftByUL & (1 | 2)))
   { // Switch; check if need to move to a new element
       case 0:
           twoBitST->elmOnUC = 0;
           return;
       case 1:
           twoBitST->elmOnUC = 1;
           return;
       case 2:
          twoBitST->elmOnUC = 2;
          return;
       case 3:
          twoBitST->elmOnUC = 3;   // On the last element
          return;
       case 4:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 0;  // Starting a new limb
          return;
       case 5:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 1; //2nd element in next limb
          return;
       case 6:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 2; //3rd element in next limb
          return;
   } // Switch; check if need to move to a new element

   return;
} // twoBitAryMoveForXElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitST to point to the previous element
\--------------------------------------------------------*/
void twoBitAryMoveBackOneElm(
  struct twoBitAry *twoBitST // array to move back in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-1 Sub-1: twoBitAryMoveBackOneElm
   '  - Moves back one element in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Move the bit above 
   switch(twoBitST->elmOnUC)
   { // Switch; check if need to move to a new element
       case 0:               // target is in previous limb
          twoBitST->elmOnUC = 3;
          --twoBitST->limbOnUCPtr; 
          return;
       case 1:
       case 2:
       case 3:
          --twoBitST->elmOnUC;
          return;
   } // Switch; check if need to move to a new element

   return;
} // twoBitAryMoveToNextElm

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBotST to point to X elements back    
\--------------------------------------------------------*/
void twoBitAryMoveBackXElm(
  struct twoBitAry *twoBitST, // To bit array to move back
  unsigned long shiftByUL  // number elements to shift back
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: twoBitAryMoveBackXElm
   '  - Moves back X elements back in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Move back by limbs
   twoBitST->limbOnUCPtr -= (shiftByUL >> 2);

   // Move the bit above 
   switch(twoBitST->elmOnUC)
   { // Switch; check if need to move to a new element
       case 0:
       // Case 0: on the first element in the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:         // Target in previous limb
                   twoBitST->elmOnUC = 3;
                   --twoBitST->limbOnUCPtr;
                   return;
               case 2:
                   twoBitST->elmOnUC = 2;
                   --twoBitST->limbOnUCPtr;
                   return;
               case 3:
                   twoBitST->elmOnUC = 1;
                   --twoBitST->limbOnUCPtr;
                   return;
           } // Switch; find out how many bits to adjust
       // Case 0: on the first element in the limb

       case 1:
       // Case 1: Was on the second element in the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:
                   twoBitST->elmOnUC = 0;
                   return;
               case 2:
                   twoBitST->elmOnUC = 3;
                   --twoBitST->limbOnUCPtr;
                   return;
               case 3:
                   twoBitST->elmOnUC = 2;
                   --twoBitST->limbOnUCPtr;
                   return;
           } // Switch; find out how many bits to adjust
       // Case 1: Was on the second element in the limb

       case 2:
       // Case 2: On the third element in the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:
                   twoBitST->elmOnUC = 1;
                   return;
               case 2:
                   twoBitST->elmOnUC = 0;
                   return;
               case 3:
                   twoBitST->elmOnUC = 3;
                   --twoBitST->limbOnUCPtr;
                   return;
           } // Switch; find out how many bits to adjust
       // Case 2: On the third element in the limb

       case 3:
       // Case 3: was on the last element of the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:
                   twoBitST->elmOnUC = 2;
                   return;
               case 2:
                   twoBitST->elmOnUC = 1;
                   return;
               case 3:
                   twoBitST->elmOnUC = 0;
                   return;
           } // Switch; find out how many bits to adjust
       // Case 3: was on the last element of the limb

       default: return; // Invalid number of elements
   } // Switch; check if need to move to a new element

   return;
} // twoBitAryMoveBackXElm

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-1 Sub-1: changeTwoBitElm
   '  - Changes a single two bit value in a two bit array.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*This adds a very small speed boost (mirco optimize)*/
    /*Get position to move bits to*/
    char valC = twoBitST->elmOnUC << 1;

    /*Set up clearing value and value to change to*/
    char clearC = 3 << valC;
    valC = newValueUC << valC;

    *twoBitST->limbOnUCPtr =
      (*twoBitST->limbOnUCPtr & (~clearC)) | valC;
    /*This clears and then sets the value to its correct
    `  postion:
    ` & (~clearC):
    `   Removes any set bits at the target location
    ` | valC
    `   sets the postion to the input value (| valC).
    */

    return;
} /*changeTwoBitElm*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    - Sets the current limb in a two-bit array to 0
\--------------------------------------------------------*/
void blankLimb(
  struct twoBitAry *twoBitST // two bit array to restart
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-1 Sub-1: blankLimb
   '  - Sets a limb in a two bit array to 0.
   '    Each limb holds four two bit elments.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   twoBitST->limbOnUCPtr = 0;
   twoBitST->elmOnUC = 0;
} // blankLimb

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Moves to the next limb in a two bit array
|  - Returns:
|    o 0: For a succesfull move
|    o 1: For memory out of bounds error
\--------------------------------------------------------*/
void moveToNextLimb(
  struct twoBitAry *twoBitST // two-bit array to move back
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-1 Sub-1: moveToNextLimb
   '  - Moves back to start of limb in a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   ++twoBitST->limbOnUCPtr;
   twoBitST->elmOnUC = 0;

   return;
} // moveToNextLimb

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o lastLimbUCPtr to point to the previous limb
\--------------------------------------------------------*/
void moveToLastLimb(
    struct twoBitAry *twoBitST // array to move back in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: Sec-1 Sub-1: moveToLastLimb
   '  - Moves to the start of the previous limb (uint8_t)
   '    in the two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   --twoBitST->limbOnUCPtr;
   twoBitST->elmOnUC = 0;

   return;
} // moveToLastLimb

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: makeTwoBitArry
   '  - Make a two bit array struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct twoBitAry *twoBitST =
     malloc(sizeof(struct twoBitAry));

   if(twoBitST == 0) return 0; // Memory error
   twoBitST->elmOnUC = 0;

   switch(blankAryBl)
   { // Switch check if making an array
     case 1:
       twoBitST->firstLimbUCPtr = 0;
       twoBitST->limbOnUCPtr = 0;
       twoBitST->lenAryUL = 0;
       return twoBitST;
   } // Switch check if making an array

   twoBitST->firstLimbUCPtr =
     calloc((lenArryUL >> 2) + 1, sizeof(uint8_t));
     // Each limb (uint8_t) has four elements, so I need
     // to divide the number of elements by 4 (>> 2)

   if(twoBitST->firstLimbUCPtr == 0) return 0;

   twoBitST->lenAryUL = ((lenArryUL >> 2) << 2) + 4;
   twoBitST->limbOnUCPtr = twoBitST->firstLimbUCPtr;

   return twoBitST;
} // makeTwoBitArray

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o dupTwoBitST to hold same pointers/values as
|      cpTwoBitST
\--------------------------------------------------------*/
void cpTwoBitPos(
  struct twoBitAry *cpTwoBitST,  // Structer to copy
  struct twoBitAry *dupTwoBitST  // Struct to copy to
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: Sec-1 Sub-1: cpTwoBitPos
   '  - copy pointers of cpTwoBitST to dupTwoBitST
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   dupTwoBitST->firstLimbUCPtr =cpTwoBitST->firstLimbUCPtr;
   dupTwoBitST->limbOnUCPtr = cpTwoBitST->limbOnUCPtr;
   dupTwoBitST->elmOnUC = cpTwoBitST->elmOnUC;
   dupTwoBitST->lenAryUL = cpTwoBitST->lenAryUL;

   return;
} // cpTwoBitPos

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: freeTwoBitAry
   '  - Frees a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(doNotFreeAryBl)
   { // switch, check if freeing the array
     case 0:
       if(stToFree->firstLimbUCPtr != 0)
         free(stToFree->firstLimbUCPtr);
   } // switch, check if freeing the array
   
   switch(twoBitOnStackBl)
   { // switch; check if freeing the twobit structer
     case 1: free(stToFree);
   } // switch; check if freeing the twobit structer

   return;
} // freeTwoBitAry

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAry to point to x elements after the start of
|      the array
\--------------------------------------------------------*/
void moveXElmFromStart(
  struct twoBitAry *twoBitST,// Array to change index for
  unsigned long shiftByUL     // Number elements to shift by
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-15 TOC: Sec-1 Sub-1: moveXElmFromStart
   '  - Change the two bit array postion to x elements
   '    from the starting position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   twoBitST->limbOnUCPtr =
      twoBitST->firstLimbUCPtr + (shiftByUL >> 2);
   twoBitST->elmOnUC = shiftByUL & (1 | 2);
   return;

   // Get whole shifts to perform
   twoBitST->limbOnUCPtr =
     twoBitST->firstLimbUCPtr + (shiftByUL >> 2);

   // Move the bit above 
   switch(shiftByUL & (1 | 2))
   { // Switch; check if need to move to a new element
       case 0:
           twoBitST->elmOnUC = 0;
           return;
       case 1:
           twoBitST->elmOnUC = 1;
           return;
       case 2:
          twoBitST->elmOnUC = 2;
          return;
       case 3:
          twoBitST->elmOnUC = 3;   // On the last element
          return;
       case 4:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 0;  // Starting a new limb
          return;
       case 5:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 1; //2nd element in next limb
          return;
       case 6:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 2; //3rd element in next limb
          return;
   } // Switch; check if need to move to a new element

   return;
} // moveXElmFromStart

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Length of the two bit array
\--------------------------------------------------------*/
unsigned long lenTwoBitAry(
  struct twoBitAry *twoBitST
    // two-bit array to get length for
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: lenTwoBitAry
   '  - Get the length of a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   return twoBitST->lenAryUL << 2;
} // lenTwoBitAry

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o The index of the elemen on in the twoBitST struct
\--------------------------------------------------------*/
unsigned long getTwoBitAryIndex(
  struct twoBitAry *twoBitST
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-17 TOC: Sec-01: getTwoBitAryIndex
   '  - Returns the index of the current two bit element
   '    on
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  return (
     ((twoBitST->limbOnUCPtr-twoBitST->firstLimbUCPtr) <<2)
    + twoBitST->elmOnUC
  ); // Return the index
} // getTwoBitAryIndex

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-18 TOC: Sec-1 Sub-1: twoBitAryMoveToNextElm
   '  - Moves to the next element in a two-bit array and
   '    returns error if out of bounds
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(twoBitST->elmOnUC == 3)
   { // If on the last element of a limb

     if(twoBitST->firstLimbUCPtr - twoBitST->limbOnUCPtr
        >= twoBitST->lenAryUL
     ) return 1; // If at end of array

   } // If on the last element of a limb
 
   // Move the bit above 
   switch(twoBitST->elmOnUC)
   { // Switch; check if need to move to a new element
       case 0:
       case 1:
       case 2:
          ++twoBitST->elmOnUC; // Move to next element
          return 0;
       case 3:            // Next element is in next limb
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 0;
          return 0;
   } // Switch; check if need to move to a new element

   return 0;
} // twoBitAryMoveToNextElm

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
  unsigned long shiftByUL  // Number elements to shift by
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-19 TOC: Sec-1 Sub-1: twoBitAryMoveForXElm
   '  - Moves forward in two bit array by x elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Get whole shifts to perform
   twoBitST->limbOnUCPtr += (shiftByUL >> 2);

   // Move the bit above 
   switch(twoBitST->elmOnUC + (shiftByUL & (1 | 2)))
   { // Switch; check if need to move to a new element
       case 0:
           twoBitST->elmOnUC = 0;
           break;
       case 1:
           twoBitST->elmOnUC = 1;
           break;
       case 2:
          twoBitST->elmOnUC = 2;
          break;
       case 3:
          twoBitST->elmOnUC = 3;   // On the last element
          break;
       case 4:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 0;  // Starting a new limb
          break;
       case 5:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 1; //2nd element in next limb
          break;
       case 6:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 2; //3rd element in next limb
          break;
   } // Switch; check if need to move to a new element

   if(twoBitST->limbOnUCPtr - twoBitST->firstLimbUCPtr >
      twoBitST->lenAryUL
   ) return 1; // Array out of bounds error

   return 0;
} // twoBitAryMoveForXElm

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-20 TOC: Sec-1 Sub-1: twoBitAryMoveBackOneElm
   '  - Moves back one element in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Move the bit above 
   switch(twoBitST->elmOnUC)
   { // Switch; check if need to move to a new element
       case 0:               // target is in previous limb
          if(twoBitST->limbOnUCPtr ==
             twoBitST->firstLimbUCPtr
          ) return 1;    // Out of bounds

          twoBitST->elmOnUC = 3;
          --twoBitST->limbOnUCPtr; 
          return 0;
       case 1:
       case 2:
       case 3:
          --twoBitST->elmOnUC;
          return 0;
   } // Switch; check if need to move to a new element

   return 2;
} // twoBitAryMoveToNextElm

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
  struct twoBitAry *twoBitST,// To bit array to move back
  unsigned long shiftByUL // number elements to shift back
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-21 TOC: Sec-1: twoBitAryMoveBackXElmBoundsCheck
   '  - Moves back X elements back in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Move back by limbs
   twoBitST->limbOnUCPtr -= (shiftByUL >> 2);

   // Move the bit above 
   switch(twoBitST->elmOnUC)
   { // Switch; check if need to move to a new element
       case 0:
       // Case 0: on the first element in the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: break; // Finished
               case 1:         // Target in previous limb
                   twoBitST->elmOnUC = 3;
                   --twoBitST->limbOnUCPtr;
                   break;
               case 2:
                   twoBitST->elmOnUC = 2;
                   --twoBitST->limbOnUCPtr;
                   break;
               case 3:
                   twoBitST->elmOnUC = 1;
                   --twoBitST->limbOnUCPtr;
                   break;
           } // Switch; find out how many bits to adjust
       // Case 0: on the first element in the limb

       case 1:
       // Case 1: Was on the second element in the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: break; // Finished
               case 1:
                   twoBitST->elmOnUC = 0;
                   break;
               case 2:
                   twoBitST->elmOnUC = 3;
                   --twoBitST->limbOnUCPtr;
                   break;
               case 3:
                   twoBitST->elmOnUC = 2;
                   --twoBitST->limbOnUCPtr;
                   break;
           } // Switch; find out how many bits to adjust
       // Case 1: Was on the second element in the limb

       case 2:
       // Case 2: On the third element in the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: break; // Finished
               case 1:
                   twoBitST->elmOnUC = 1;
                   break;
               case 2:
                   twoBitST->elmOnUC = 0;
                   break;
               case 3:
                   twoBitST->elmOnUC = 3;
                   --twoBitST->limbOnUCPtr;
                   break;
           } // Switch; find out how many bits to adjust
       // Case 2: On the third element in the limb

       case 3:
       // Case 3: was on the last element of the limb
           switch(shiftByUL & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: break; // Finished
               case 1:
                   twoBitST->elmOnUC = 2;
                   break;
               case 2:
                   twoBitST->elmOnUC = 1;
                   break;
               case 3:
                   twoBitST->elmOnUC = 0;
                   break;
           } // Switch; find out how many bits to adjust
       // Case 3: was on the last element of the limb

       default: return 2; // Invalid number of elements
   } // Switch; check if need to move to a new element

   if(twoBitST->limbOnUCPtr < twoBitST->firstLimbUCPtr)
     return 1; // Array out of bounds error

   return 0;
} // twoBitAryMoveBackXElm

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-22 TOC: Sec-1 Sub-1: moveToNextLimb
   '  - Moves back to start of limb in a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   ++twoBitST->limbOnUCPtr;
   twoBitST->elmOnUC = 0;

   if(twoBitST->limbOnUCPtr - twoBitST->firstLimbUCPtr 
      > twoBitST->lenAryUL
   ) return 1;

   return 0;
} // moveToNextLimb

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-23 TOC: Sec-1 Sub-1: moveToLastLimb
   '  - Moves to the start of the previous limb (uint8_t)
   '    in the two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(twoBitST->limbOnUCPtr == twoBitST->firstLimbUCPtr)
   { // If I am on the first limb
     // Check If I can shift one to the first element
     if(twoBitST->elmOnUC == 0) return 2;

     twoBitST->elmOnUC = 0;
     return 1;
   } // If I am on the first limb

   --twoBitST->limbOnUCPtr;
   twoBitST->elmOnUC = 0;

   return 0;
} // moveToLastLimb

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
  unsigned long shiftByUL // Number elements to shift by
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-24 TOC: Sec-1 Sub-1: moveXElmFromStart
   '  - Change the two bit array postion to x elements
   '    from the starting position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Get whole shifts to perform
   twoBitST->limbOnUCPtr =
     twoBitST->firstLimbUCPtr + (shiftByUL >> 2);

   // Move the bit above 
   switch(shiftByUL & (1 | 2))
   { // Switch; check if need to move to a new element
       case 0:
           twoBitST->elmOnUC = 0;
           break;
       case 1:
           twoBitST->elmOnUC = 1;
           break;
       case 2:
          twoBitST->elmOnUC = 2;
          break;
       case 3:
          twoBitST->elmOnUC = 3;   // On the last element
          break;
       case 4:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 0;  // Starting a new limb
          break;
       case 5:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 1; //2nd element in next limb
          break;
       case 6:
          ++twoBitST->limbOnUCPtr;
          twoBitST->elmOnUC = 2; //3rd element in next limb
          break;
   } // Switch; check if need to move to a new element

   if(twoBitST->limbOnUCPtr - twoBitST->firstLimbUCPtr >
      twoBitST->lenAryUL
   ) return 1; // Array out of bounds error

   return 0;
} // moveXElmFromStart

