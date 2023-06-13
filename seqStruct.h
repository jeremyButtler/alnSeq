/*#########################################################
# Name: seqStruct
# Use:
#  - Holds the seqStruct and functions to manipulate
#    the seqStruct strucuter
# Libraries:
# C Standard Libraries:
#  - <stdint.h>
#  - <stdlib.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o struct-01: seqStruct
'    - Holds sequence and length of a input sequence
'  o fun-01 initSeqST:
'     - Sets vlues in seqST to zero
'  o fun-02 setSeqSTSequence:
'     - Sets vlues in seqST to zero
'  o fun-03 setSeqSTQScore:
'    - Sets the Q-score for a sequence and finds the
'  o fun-04 addStartEndToSeqST:
'     - Sets the start and ending corrdinates of a region
'       of interest in a sequence
'  o fun-05 addSeqId:
'    - Adds the sequence id to the seqST struct
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef SEQSTRUCT_H
#define SEQSTRUCT_H

#include <stdint.h>
#include <stdlib.h>

/*--------------------------------------------------------\
| Struct-01: seqStruct
|  - Holds sequence and length of a input sequence
\--------------------------------------------------------*/
typedef struct seqStruct
{ // refQueryStruct
  char *idCStr;          // Id of th sequence
  uint32_t lenIdUC;      // Length of the sequence id
  uint32_t lenIdBuffUI;  // Lenght of Id buffer

  char *seqCStr;          // Sequence
  uint32_t lenSeqUI;      // Length of the sequence
  uint32_t lenSeqBuffUI;  // Lenght of sequence buffer

  char *qCStr;           // q-score entry
  uint32_t lenQUI;       // Length of the Q-score
  uint32_t lenQBuffUI;   // Length of Q-score buffer

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
|    o qST to point to qCStr and hold the length of qCStr
|  - Frees
|    o The old q-score (if there was an old q-score) if
|      freeOldQBl is not  0
| Note:
|  - This is a shallow copy (pointer only), so do not free
|    qCStr
|  - Make sure you have a pointer to seqST->qCStr if
|    you did set freeOldQBl to 0. Otherwise their is
|    risk of lossing a handle on a pointer.
\--------------------------------------------------------*/
void setSeqSTQScore(
  char *qCStr,          // Sequence to add to structure
  char freeOldQBl,      // 0: do not free old sequence
  uint32_t lenQUI,      // New q entry length (0 to find)
  struct seqStruct *seqST // Struct to add sequence to
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-01: setSeqSTQScore
   '  - Sets the Q-score for a sequence and finds the
   '    length if 0 is input
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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-01: addSeqId
   '  - Adds the sequence id to the seqST struct
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
