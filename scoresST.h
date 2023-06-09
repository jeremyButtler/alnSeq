/*#########################################################
# Name: scoreST
# Use:
#  - Hols the score structure,which holds the index and
#    score of a single cell in a Needleman-Wunsch or
#    Smith-Waterman direction matrix
# Libraries:
# C Standard Libraries:
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o struct-01: scoresStruct
'    - Holds the score for a single direction matrix
'      positon
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

#ifndef SCOREST_H
#define SCORESST_H

/*--------------------------------------------------------\
| Struct-01: scoresStruct
|  - Holds the score for a single direction matrix positon
\--------------------------------------------------------*/
typedef struct scoresStruct
{ // scoresStruct
  unsigned long indexUL;// Index in dirtion matrix of score
  long scoreL;          // Score of positoin
}scoresStruct;

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: sortScores
   '  - Sorts an array of scores structs by score using
   '    shell short.  Order is greatest to least.
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '  o fun-01 sec-02:
   '    - Build the sub array sizes for each round of shell
   '      sort (sub array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o firstScoreST to hold the values from secScoreST and
|      secScoreST to hold values from firstScoreST
\--------------------------------------------------------*/
void swapScoreSTs(
  struct scoresStruct *firstScoreST, // Has values to swap
  struct scoresStruct *secScoreST    // Has values to swap
);  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-02 TOC: Sec-01 Sub-01: swapScoreSTs
    '  - Swaps values in two score structures
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o scoreST to have values set to 0
\--------------------------------------------------------*/
void initScoresST(
  struct scoresStruct *scoreST
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-0 Sub-0: initScoresST
   '  - Sets scores in a scores struture to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o Frees the score struct if requested
\--------------------------------------------------------*/
void freeScoresST(
  struct scoresStruct *scoreST, // struct to free
  char stackBl  // 1: structure is on the stack
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-0 Sub-0: freeScoresST
   '  - Frees score structure if on heap, else does nothing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreST to have the new index
\--------------------------------------------------------*/
void changeScoreSTIndex(
  struct scoresStruct *scoreST, // has index to change
  unsigned long newIndexUL      // index to swap to
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-0 Sub-0: changeScoreSTIndex
   '  - Changes the index value in a scores structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o scoreST to have the new score
\--------------------------------------------------------*/
void changeScoreSTScore(
  struct scoresStruct *scoreST, // has index to change
  long newScoreL               // score to swap to
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-0 Sub-0: changeScoreSTScore
   '  - Changes the score value in a scores structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
