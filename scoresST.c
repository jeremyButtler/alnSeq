/*#########################################################
# Name: scoreST
# Use:
#  - Hols the score structure,which holds the index and
#    score of a single cell in a Needleman-Wunsch or
#    Smith-Waterman direction matrix
# Libraries:
# C Standard Libraries:
#########################################################*/

#include "sortScores.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 sortScores:
'    - Sorts an array of scores structs by score using
'      shell short.  Order is greatest to least.
'  o fun-02 swapScoreSTs:
'    - Swaps values in two score structures
'  o fun-03 initScoresST:
'    - Sets scores in a scores struture to 0
'  o fun-04 freeScoresST:
'    - Frees score structure if on heap, else does nothing
'  o fun-05 changeScoreSTIndex:
'    - Changes the index value in a scores structure
'  o fun-06 changeScoreSTScore:
'    - Changes the score value in a scores structure
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Outptu:
|  - Modifies
|    o sorts the scoresST array between the index of the
|      first element to the index of the last element
\--------------------------------------------------------*/
void sortScores(
   struct scoresStruct **scoresST,//array of socres to sort
   unsigned long firstElmUL,//Index of 1st element to sort
   unsigned long lastElmUL, //Index of last element to sort
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: sortScores
   '  - Sorts an array of scores structs by score using
   '    shell short. Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-01 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  // Number of elements to sort
  unsigned long numElmUL = lastElmUL - firstElmUL;

  // Number of sorting rounds
  unsigned long rndsUL = (numElmUL / 3);
  unsigned long subAryUL[rndsUL];
    // Array of sizes for the sub arrays in shell sort

  unsigned long lastScoreUL = 0;
  unsigned long tmpScoreUL = 0;
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-01 Sec-02:
  ^  - Build the sub array sizes for each round of shell
  ^     sort (sub array)
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  // Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] + 1
  subAryUL[0] = 1; // Initialzie first array

  for(unsigned long cntUL = 1; cntUL < rndsUL; ++cntUL)
    subAryUL[cntUL] = (3 * subAryUL[cntUL - 1]) + 1;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-01 Sec-03:
  ^  - Sort the scores structure array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  while(rndsUL > 0)
  { // loop trhough all sub arrays sort the subarrays
    --rndsUL; // Move to the next subarray round
      // Account for index 1

    for(
      unsigned long elmUL = 0;
      elmUL < subAryUL[rndsUL];
      ++elmUL
    ) { // For each element in the subarray

      for(
        unsigned long subAryOnUL = firstElmUL + elmUL;
        subAryOnUL + subAryUL[rndsUL] <= lastElmUL;
        subAryOnUL += subAryUL[rndsUL];
      ) { // Loop; swap each nth element of the subarray

        if(*(scoresST + subAryOnUL)->scoreL <
             *(scoreST+subAryOnUL+subAryUL[rndsUL])->scoreL
        ) { // If I need to swap elements

          swapScoreSTs(
            *(scoresST + subAryOnUL),
            *(scoreST + subAryOnUL + subAryUL[rndsUL])
          ); // Swap scores around

          lastScoreUL = subAryOnUL;
          tmpScoreUL = lastScoreUL - subAryUL[rndsUL];

          while(
            tmpScoreUL >= firstElmUL &&
            *(scoresST + lastScoreUL)->scoreL >
              *(scoreST + tmpScoreUL)->scoreL
          ) { // loop; move swapped element back
            swapScoreSTs(
              *(scoresST + lastScoreUL),
              *(scoreST + lastScoreUL - subAryUL[rndsUL])
            ); // Swap socres around

            tmpScoreUL -= subAryUL[rndsUL];
          } // loop; move swapped element back
        } // If I need to swap elements
      } // Loop; swap each nth element of the subarray
    } // For each element in the subarray
  } // loop trhough all sub arrays sort the subarrays

  return; // Finshed sorting the array
} // sortScores

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o firstScoreST to hold the values from secScoreST and
|      secScoreST to hold values from firstScoreST
\--------------------------------------------------------*/
void swapScoreSTs(
  struct scoresStruct *firstScoreST, // Has values to swap
  struct scoresStruct *secScoreST    // Has values to swap
){  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-02 TOC: Sec-01 Sub-01: swapScoreSTs
    '  - Swaps values in two score structures
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  unsigned long indexUL = firstScore->indexUL;
  long scoreL = firstScore->scoreL;

  firstScoreST->indexUL = secScoreST->indexUL;
  firstScoreST->scoreL = secScoreST->scoreL;

  secScoreST->indexUL = indexUL;
  secScoreST->scoreL = scoreL;
} // swapScoreSTs

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o scoreST to have values set to 0
\--------------------------------------------------------*/
void initScoresST(
  struct scoresStruct *scoreST // Struct to initalize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-0 Sub-0: initScoresST
   '  - Sets scores in a scores struture to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  scoreST->indexUL = 0;
  scoreST->scoreL = 0;
  return;
} // initScoresST

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o Frees the score struct if requested
\--------------------------------------------------------*/
void freeScoresST(
  struct scoresStruct *scoreST, // struct to free
  char stackBl  // 1: structure is on the stack
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-0 Sub-0: freeScoresST
   '  - Frees score structure if on heap, else does nothing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  if(stackBl & 1) free(scoreST);
  return;
} // freeScoresST

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreST to have the new index
\--------------------------------------------------------*/
void changeScoreSTIndex(
  struct scoresStruct *scoreST, // has index to change
  unsigned long newIndexUL      // index to swap to
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-0 Sub-0: changeScoreSTIndex
   '  - Changes the index value in a scores structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  scoreST->indexUL = newIndexUL;
  return;
} // changeScoreSTIndex

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreST to have the new score
\--------------------------------------------------------*/
void changeScoreSTScore(
  struct scoresStruct *scoreST, // has index to change
  long newScoreL               // score to swap to
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-0 Sub-0: changeScoreSTScore
   '  - Changes the score value in a scores structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  scoreST->scoreL = newScoreL;
  return;
} // changeScoreSTScore
