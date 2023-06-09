/*#########################################################
# Name: seqStruct
# Use:
#  - Holds the seqStruct and functions to manipulate
#    the seqStruct strucuter
# Libraries:
# C Standard Libraries:
#  - <stdint.h>
#########################################################*/

#include "seqStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
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

/*--------------------------------------------------------\
| Output:
|  - Modifies
|    o all values in seqST to be 0
\--------------------------------------------------------*/
void initSeqST(
  struct seqStruct *seqST // Struct to initialize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-01: initSeqST
   '  - Sets vlues in seqST to zero
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   seqST->idCStr = 0;
   seqST->lenIdUC = 0;
   seqST->seqCStr = 0;
   seqST->lenSeqUI = 0;
   seqST->offsetUI = 0;
   seqST->endAlnUI = 0;

   return;
} // initSeqST

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-01: setSeqSTSequence
   '  - Sets vlues in seqST to zero
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(seqST->seqCStr != 0 && freeOldSeqBl != 0)
     free(seqST->seqCStr);

   seqST->seqCStr = seqCStr;

   if(lenSeqUI > 0) seqST->lenSeqUI = lenSeqUI;

   else
   { // Else need to find the lenght of the sequence
     seqST->lenSeqUI = 0;

     while(*seqCStr != '\0') 
     { // While not at the end of the sequence
       ++seqCStr;
       ++seqST->lenSeqUI;
     } // While not at the end of the sequence
   } // Else need to find the lenght of the sequence

   return;
} // setSeqSTSequence

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-01: addStartEndToSeqST
   '  - Sets the start and ending corrdinates of a region
   '    of interest in a sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   seqST->offsetUI = startTargetUI;
   seqST->endAlnUI = endTargetUI;

   return;
} // addStartEndToSeqST

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

   if(seqST->idCStr != 0 && freeOldIdBl != 0)
     free(seqST->idCStr);

   seqST->idCStr = idCStr;
   
   if(lenIdUC == 0)
   { // If I need ot find the length of the infput id
     seqST->lenIdUC = 0;

     while(*idCStr != 0)
     { // While I have the id lenght to find
       ++idCStr;
       ++seqST->lenIdUC;
     } // While I have the id lenght to find
   } // If I need ot find the length of the infput id

   else seqST->lenIdUC = lenIdUC;

   return;
} // addSeqId

/*--------------------------------------------------------\
| Output:
|  - Frees
|    o seqST from memory
| Notes:
|  - You will have to set the pointer to seqST to 0
|  - Make sure you have a pointer to seqST->seqCStr if
|    you did set freeSeqBl to 0. Otherwise you will loose
|    your handle to the data in seqCStr
\--------------------------------------------------------*/
void freeSeqST(
  struct seqStruct *seqST, // Struct to free
  char freeSeqBl,  // 0: do not free seqCStr
  char freeIdBl,   // 0: do not free the id
  char heapBl     // 0: seqST on stack only free variables
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-01: freeSeqST
   '  - Frees the seqST strucuter
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(freeSeqBl != 0 && seqST->seqCStr != 0)
     free(seqST->seqCStr);
   
   if(freeIdbl != 0 && seqST->idCStr != 0)
     free(seqST->idCStr);

   seqST->seqCStr = 0;
   seqST->idCStr = 0;

   if(heapBl != 0) free(seqST);

   return;
} // freeSeqST
