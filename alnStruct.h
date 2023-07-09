/*#########################################################
# Name: alnStruct
# Use:
#  - Holds the alingment structure that stores the
#    alignment. This also includes the functions needed
#    to maintain and print out the alignment
# Libraries:
#  - "seqStruct.h"
#  - "scoresST.h"
#  - "alnSetStruct.h"
#  - "generalAlnFun.h"
#  o "twoBitArrays.h"
#  o "alnSeqDefaults.h"
# C Standard Libraries:
#  o <stdio.h>
#  o <stdlib.h>
#  o <stdint.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o struct-01: alnStruct
'    - Holds the alignment array, which tells if a match,
'      SNP, del, ins, or masking was at a postion
'  o fun-01 cnvtAlnAryToSeq:
'    - Uses an array of error types and a c-string with
'      the a sequence to make one part of an alignment
'  o fun-02 alnAryToLetter:
'    - Converts an alignment error array from my alignment
'      algorithms into an array of letters
'      (I = insertion, D = deletion, = = match, X = snp)
'  o fun-03 cnvtDirMatrixToAlnAry:
'    - Builds an anlignment array for the input direction
'  o fun-04 printAln:
'    - Prints out an alignment
'  o fun-05 initAlnST:
'    - Initalize all values in alnST to 0
'  o fun-06 freeAlnST:
'    - Frees alnST and all variables in alnST
'  o fun-07 addAlnSTArray:
'    - Adds an alignment array to an alnStruct structure
'  o fun-08 lnSTAddNewCode:
'    - Convert a two bit element into an alignment code.
'      Matches are checked in fun-02 (alnAryToLetter)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef ALNSTRUCT_H
#define ALNSTRUCT_H

#include "seqStruct.h" // Need to set up
#include "scoresST.h"
#include "generalAlnFun.h"

// Flags for the alignment array (struct-01)
#define defDelFlag 1    // deletion
#define defInsFlag 2    // insertion
#define defBaseFlag 4   // match or snp
#define defSoftQueryFlag 8 // Softmask a query base
#define defSoftRefFlag 16  // Softmask a reference base

#define defAlnEndFlag 0
#define defGapFlag 1
#define defSnpFlag 2
#define defMatchFlag 3

/*--------------------------------------------------------\
| Struct-01: alnStruct
|  - Holds the alignment array, which tells if a match,
|    SNP, del, ins, or masking was at a postion
\--------------------------------------------------------*/
typedef struct alnStruct
{ // alnStruct
  uint8_t *alnAryUC;    // Array of alignment types
                        // See definitions above for list
  uint32_t lenAlnAryUI; // Length of alnAryUC
  uint32_t numBasesUI;  // Number of bases in alignment
  uint32_t numAlnBasesUI; // Number of aligned bases

  // Starting and ending coordiantes for the alignment
  uint32_t refStartUI;
  uint32_t refEndUI;
  uint32_t queryStartUI;
  uint32_t queryEndUI;
}alnStruct;

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o Heap alloacted C-string with alignment for the
|      input sequence
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
char * cnvtAlnAryToSeq(
    struct seqStruct *seqST, // Has sequence to work with
    char queryBl,            // 1: working on query; 0: ref
    struct alnStruct *alnST  // Has alignment array
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: cnvtAlnQueryAryToSeq
   '  - Uses an array of error types and a c-string with
   '    the a sequence to make one part of an alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Algiment array in alnST to have letters instead of
|      codes
| Note:
|  - Only call this function after you are done with alnST
\--------------------------------------------------------*/
void alnAryToLetter(
    struct seqStruct *refST,  // Ref sequence for detecting
                              // matches; Input 0 to ignore
    struct seqStruct *queryST,// query seq for detecting
                              // matches Input 0 to ignore
    struct alnStruct *alnST // Has alignment array to
                            // convert to human readable
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: cnvtAlnErrAryToLetter
   '  - Converts an alignment array from my alignment
   '    algorithms into an array of letters
   '    (I = insertion, D = deletion, = = match, X = snp)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o alnStruct with the alingment array
|    o 0 if had memory allocation error
\--------------------------------------------------------*/
struct alnStruct * cnvtDirMatrixToAlnAry(
    struct seqStruct *refST,  // Has reference seq & length
    struct seqStruct *queryST,   // Query sequence & length
    struct twoBitAry *dirMatrxST, // Direction matrix
    struct scoresStruct *scoreST, //has index of best score
    char softMaskBl       // 1: Apply soft masking
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: cnvtDirMatrixToAlnAry
   '  - Builds an anlignment array for the input direction
   '    matrix
   '  o fun-03 sec-01:
   '    - VariAble declerations
   '  o fun-03 sec-02:
   '    - BuilD the alignment array
   '  o fun-03 sec-03:
   '    - Clean up and add softmasking to start
   '  o fun-03 sec-04:
   '    - InveRt the error type array
   '  o fun-03 sec-05:
   '    - Add softmasking to the end and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints out the input alingment
\--------------------------------------------------------*/
void printAln(
  FILE *outFILE,      // File to print alingment to
  char *queryIdCStr,  // Id/name of query sequence
  char *queryAlnCStr, // Alinged query sequence
  char *refIdCStr,    // Id/name of reference sequence
  char *refAlnCStr,   // Aligned reference sequence
  long scoreL,        // Score of the alignment
  unsigned short lineWrapUS,
    // Number of characters per line (minimum is 42)
  struct alnStruct *alnST // alignment array in human form
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: printAln
   '  - Prints out the alignment
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Print out the header  for the alignment
   '  o fun-04 sec-04:
   '    - Print out tail of the alingment (missed by loop)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o All variables in alnST to be 0
| Note:
|  - This does not free alnAryUc, so only call this for
|    new alnST structures
\--------------------------------------------------------*/
void initAlnST(
  struct alnStruct *alnST // Strucutre to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-01: initAlnST
   '  - Initalize all values in alnST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o alnST, including all variables in in alnST
|    0 Or if heapBl is 0; frees all variables in alnST, but 
|      does not free alnST
\--------------------------------------------------------*/
void freeAlnST(
  struct alnStruct *alnST, // Strucutre to free
  char heapBl // 0: free variables in alnST, but keep alnST
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-01: freeAlnST
   '  - Frees alnST and all variables in alnST
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o alnST->alnAryUC to have an alignment array
|    o alnST->lenAlnAryUI to have the length of the
|      alignment array
|  - Returns
|    o 0 for no errors
|    o 64 for memory allocation errors
\--------------------------------------------------------*/
unsigned char addAlnSTArray(
  struct alnStruct *alnST, // add new array to
  unsigned long  lenAryUL  // Size of new array
    // This should be reference length + query length
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: addAlnSTArray
   '  - Adds an alignment array to an alnStruct structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o alnST to hold the alignment direction and incurment
|      both of its counters
\--------------------------------------------------------*/
void alnSTAddNewCode(
  struct alnStruct *alnST, // Has alignment array to add to
  uint8_t twoBitElmUC      // Alignment code to add 
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: alnSTAddNewCode
   '  - Convert a two bit element into an alignment code.
   '    Matches are checked in fun-02 (alnAryToLetter)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
