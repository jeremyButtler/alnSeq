/*#########################################################
# Name: seqStruct
# Use:
#  - Holds the seqStruct and functions to manipulate
#    the seqStruct strucuter
# Libraries:
# C Standard Libraries:
#  - <stdint.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o struct-01: seqStruct
'    - Holds sequence and length of a input sequence
'  o fun-01 initSeqST:
'     - Sets vlues in seqST to zero
'  o fun-02 setSeqSTSequence:
'     - Sets vlues in seqST to zero
'  o fun-03 addStartEndToSeqST:
'     - Sets the start and ending corrdinates of a region
'       of interest in a sequence
'  o fun-04 addSeqId:
'    - Adds the sequence id to the seqST struct
'  o fun-05 freeSeqST:
'     - Frees the seqST strucuter
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef SEQSTRUCT_H
#define SEQSTRUCT_H

#include <stdint.h>

/*--------------------------------------------------------\
| Struct-01: seqStruct
|  - Holds sequence and length of a input sequence
\--------------------------------------------------------*/
typedef struct seqStruct
{ // refQueryStruct
  char *idCStr;          // Id of th sequence
  unsigned char lenIdUC; // Length of the sequence id
  char *seqCStr;         // Sequence
  uint32_t lenSeqUI;     // Length of the sequence
  uint32_t offsetUI;     // Offset for an alignment
  uint32_t endAlnUI;     // Marks end of alignment
}seqStruct;

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o all values in seqST to be 0
\--------------------------------------------------------*/
void initSeqST(
  struct seqStruct *seqST // Struct to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-01: initSeqST
   '  - Sets vlues in seqST to zero
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o seqST to point to seqCStr and hold the length of
|      seqCStr
|  - Frees
|    o The old sequence (if there was an old sequence) if
|      freeOldSeqBl is not  0
| Note:
|  - This is a shallow copy (pointer only), so do not free
|    seqCStr
|  - Make sure you have a pointer to seqST->seqCStr if
|    you did set freeOldSeqBl to 0. Otherwise their is
|    risk of lossing a handle on a pointer.
\--------------------------------------------------------*/
void setSeqSTSequence(
  char *seqCStr,          // Sequence to add to structure
  uint32_t lenSeqUI,      // Length of seq (0 if not known)
  char freeOldSeqBl,      // 0: do not free old sequence
  struct seqStruct *seqST // Struct to add sequence to
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-01: setSeqSTSequence
   '  - Sets vlues in seqST to zero
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o seqST to have the start and end corrdiantes of the
|      target region in the sequence
\--------------------------------------------------------*/
void addStartEndToSeqST(
  uint32_t startTargetUI, // Start of region of intreset
  uint32_t endTargetUI,   // End of region of interest
  struct seqStruct *seqST // Struct to add corrdinates to
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-01: addStartEndToSeqST
   '  - Sets the start and ending corrdinates of a region
   '    of interest in a sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o seqST to hold the input id
| Note:
|  - This is a shallow copy (pointer only), so do not free
|    idCStr
|  - Make sure you have a pointer to the old seqST->idCStr
|    if you set freeOldIdBl to 0. Otherwise you will lose
|    your handle on the old id.
\--------------------------------------------------------*/
void addSeqId(
  char *idCStr,           // Sequence id to add in
  unsigned char lenIdUC,  // Length of id (0 if unknown)
  char freeOldIdBl,       // 0: Do not free old seq id
  struct seqStruct *seqST // Struct to initialize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-01: addSeqId
   '  - Adds the sequence id to the seqST struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o seqST from memory
| Notes:
|  - You will have to set the pointer to seqST to 0
|  - Make sure you have a pointer to seqST->seqCStr, and
|    seqST->idCStr if you did set freeSeqBl or freeIdbl to
|    0. Otherwise you will loose your handle to the data
\--------------------------------------------------------*/
void freeSeqST(
  struct seqStruct *seqST, // Struct to free
  char freeSeqBl,  // 0: do not free seqCStr
  char freeIdBl,   // 0: do not free the id
  char heapBl     // 0: seqST on stack only free variables
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-01: freeSeqST
   '  - Frees the seqST strucuter
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
