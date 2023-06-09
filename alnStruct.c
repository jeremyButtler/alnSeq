/*#########################################################
# Name: alnStruct
# Use:
#  - Holds the alingment structure that stores the
#    alignment. This also includes the functions needed
#    to maintain and print out the alignment
# Libraries:
#  - "twoBitArrays.h"
#  - "seqStruct.h"
#  - "scoresST.h"
#  - "alnSetStruct.h"
#  o "alnSeqDefaults.h"
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdio.h>
#  o <stdint.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
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
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define "alnStruct.h"

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Heap alloacted C-string with alignment for the
|      input sequence
|    o 0 for memory allocation errors
\--------------------------------------------------------*/
char * cnvtAlnAryToSeq(
    struct seqSturct *seqST, // Has sequence to work with
    char queryBl,            // 1: working on query; 0: ref
    struct alnStruct *alnST  // Has alignment array
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: cnvtAlnQueryAryToSeq
   '  - Uses an array of error types and a c-string with
   '    the a sequence to make one part of an alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *baseCStr = seqST->seqCStr + seqST->offsetUI - 1;
   char *tmpBaseCStr = 0;
   char *seqAlnCStr = 0;
   uint8_t *errUCPtr = alnST->alnAryUC;

   seqAlnCStr = malloc(sizeof(char)* alnST->lenErrAryUI+2);

   if(seqAlnCStr == 0) return 0;  // memory error

   tmpBaseCStr = seqAlnCStr;

   while(*errUCPtr != 0)
   { // While I have bases or an alignment to copy over
       switch(*errUCPtr)
       { // Switch; check the error type
           case defDelFlag:                    // deletion
               if(queryBl & 1) *tmpBaseCStr = '-';
               else
               { // Else dealing with a reference sequence
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } // Else dealing with a reference sequence

               break;

           case defInsFlag:                   // insertion
               if(!(queryBl & 1)) *tmpBaseCStr = '-';

               else
               { // Else dealing with a query sequence
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } // Else dealing with a query sequence

               break;

           case defSoftQueryFlag:
               if(queryBl != 0)
               { // If softmasking a query region
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } // If softmasking a query region

               // Else working on the reference sequence
               else *tmpBaseCStr = '-';
               break;


           case defSoftRefFlag:    // Soft masked ref base
               if(queryBl == 0)
               { //Else if is reference base; soft masking
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } //Else if is reference base; soft masking

               // Else working on a query sequence
               else *tmpBaseCStr = '-';
               break;

           // Is bot a query and reference soft mask
           case defSoftQueryFlag + defSoftRefFlag:
               *tmpBaseCStr = *baseCStr;
               ++baseCStr;
               break;

           case defBaseFlag:                // (match/snp)
               *tmpBaseCStr = *baseCStr;
               ++baseCStr;
               break;
       } // Switch; check the error type

       // Move to the next error type
       ++tmpBaseCStr;
       ++errUCPtr;
   } // While I have bases or an alignment to copy over

  *tmpBaseCStr = '\0'; // Make into a c-string
  return seqAlnCStr;
} // cnvtAlnQueryAryToSeq

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o Algiment array in alnST to have letters instead of
|      codes
| Note:
|  - Only call this function after you are done with alnST
\--------------------------------------------------------*/
void cnvtAlnErrAryToLetter(
    char *refSeqCStr,       // Ref sequence for detecting
                            //  matches; Input 0 to ignore
    char *querySeqCStr,     // query sequence for detecting
                            // matches Input 0 to ignore
    struct alnStruct *alnST // Has alignment array to
                            // convert to human readable
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: cnvtAlnErrAryToLetter
   '  - Converts an alignment array from my alignment
   '    algorithms into an array of letters
   '    (I = insertion, D = deletion, = = match, X = snp)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    uint8_t alnAryUC = alnST->alnAryUC;

    while(*alnAryUC != 0)
    { // While have codes to convert
        switch(*alnAryUC)
        { // Switch; check the error type
            case defDelFlag:            // deletion
                *alnAryUC = 'D';
                ++refSeqCStr;    // Move to next ref base
                break;

            case defInsFlag:            // insertion
                *alnAryUC = 'I';
                ++querySeqCStr; // Move to next query base
                break;

            case defSoftQueryFlag: // Query base sofmasked
                *alnAryUC = 's';// mark as query masked
                ++querySeqCStr; // Move to next query base
                break;

            case defSoftRefFlag: // Ref base softmasked
                *alnAryUC = 'P';// Mark as ref softmask
                ++refSeqCStr;   // Move to next query base
                break;

            case defSoftRefFlag + defSoftQueryFlag: 
                *alnAryUC = 'S'; // Mark as both masked
                ++querySeqCStr; // Move to next query base
                ++refSeqCStr;   // Move to next ref base
                break;

            case defBaseFlag:             // Match or snp
                // Check if am calling matches/snps
                if(refSeqCStr == 0 || querySeqCStr == 0)
                    *alnAryUC = 'X';

                else if(
                  checkIfBasesMatch(
                      querySeqCStr,
                      refSeqCStr
                    ) == 0) *alnAryUC = 'X';    // snp

                else *alnAryUC = '=';   // match

                ++refSeqCStr;  // Move to next ref base
                ++querySeqCStr; // Move to next query base
                break;
        } // Switch; check the error type

        ++alnAryUC;
    } // While have codes to convert

    return;
} // cnvtAlnErryAryToLetter

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: getAlnAry
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

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-01:
  ^  - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  char *bestQueryCStr =
      queryST->seqCStr
    + (   scoreST->indexUL
        / (queryST->lenRefUI - queryST->offsetUI)
      )
    - 1;
    // bestScoreUI / lenRefUI gives the number of rows
    // till the best score. -1 accounts for the  the
    // empty row
  char *bestRefCStr =
      refST->seqCStr
    + (   scoreST->indexUL
        % (refST->lenRefUI - refST->offsetUI)
      )
    - 1;
    // bestScoreUI % lenRefUI gives the number of columns
    // I am in on the row of the best score. -1 accounts
    // for the empty column

  // These recored the ending bases
  char *endQueryAlnCStr = bestQueryCStr;
  char *endRefAlnCStr = bestRefCStr;

  char *tmpQueryCStr = 0;
  char *tmpRefCStr = 0;

  struct alnStruct *alnST = 0;

  // The alignment arrary starts out at the end and needs
  // to be reorianted (inverted). These variables are for
  // this
  uint8_t *startUCPtr = 0;
  uint8_t *endUCPtr = 0;
  uint8_t swapUC = 0;

  alnST = malloc(sizeof(alnST));

  if(alnST == 0) return 0;

  alnST->lenAlnAryUI = lenQueryUI + lenRefUI;
  alnST->alnAryUC = calloc(*lenErrAryUI, sizeof(uint8_t));
    // calloc initializes all values to 0, so I do not
    // need to mark the end of the alignment array

  if(alnST->alnAryUC == 0)
  { // If had a memory error
     freeAlnST(alnST);
     return 0; // Memory error
  { // If had a memory error

  // Get the ending position for the alignment
  alnST->refEndUI = endRefAlnCStr - refCStr;
  alnST->queryEndUI = endQueryAlnCStr - queryCStr;

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-02:
  ^  - Build the alignment array
  ^  o fun-03 sec-5 sub-1:
  ^    - Find the best score cell
  ^  o fun-03 sec-5 sub-2:
  ^    - Find the best path
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*******************************************************\
  * Fun-03 Sec-5 Sub-1:
  *  - Find the cell with the best score
  \*******************************************************/

  // move to best score using the best scores index
  moveXElmFromStart(dirMatrxST, scoreST->indexUL);

  // get direction for the best element
  bitElmUC = getTwoBitAryElm(dirMatrxST);

  /*******************************************************\
  * Fun-03 Sec-5 Sub-2:
  *  - Find the best path
  \*******************************************************/

  while(bitElmUC != defMoveStop)
  { // While I have more bases in the alignment
    lastBitElmUC = bitElmUC; // Keep track of last bit

    switch(bitElmUC)
    { // Switch: check what the next base is in the sequence
      case defMoveUp:
        twoBitAryMoveBackXElm(dirMatrxST, lenRefUL + 1);
          // lenRefUL (index 1) cells per row; need + 1 to
          // get to the cell above
        *(alnErrAryUC + numErrUI) = defInsFlag;
        --bestQueryCStr; // insertion only query has a base
        break;

      case defMoveDiagnol:
        *(alnErrAryUC + numErrUI) = defBaseFlag;
        twoBitAryMoveBackXElm(dirMatrxST, lenRefUL + 2);
            // lenRefUL (index 1) cells per row; need + 2
            // to get to the next diagnol cell
        --bestQueryCStr;// match/snp query & ref have bases
        --bestRefCStr;
        break;

      case defMoveLeft:
        twoBitAryMoveBackOneElm(dirMatrxST);
        *(alnErrAryUC + numErrUI) = defDelFlag; 
        --bestRefCStr; // deletion; reference only has base
      break;
    } // Switch: check what the next base is in the sequence

      // Move back to the selected base

      // Get the next direction to move
      bitElmUC = getTwoBitAryElm(dirMatrixST);

      ++alnST->numAlnBasesUI; // Add up all aligned bases
  } // While I have more bases to add to the path

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-03:
  ^  - Clean up and add softmasking to start
  ^  o fun-03 sec-06 sub-01:
  ^    - clean up
  ^  o fun-03 sec-06 sub-02: 
  ^    - Add soft masking for the starting bases
  ^  o fun-03 sec-06 sub-03:
  ^    - Mark end of alignment & clean up indels at end
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  /*******************************************************\
  * Fun-03 Sec-06 Sub-01:
  *  - Find starting position and clean up
  \*******************************************************/

  alnST->numBasesUI = alnST->numAlnBasesUI;

  alnST->refStartUI =
      alnST->offsetUI
    + (endRefAlnCStr - bestRefCStr);

  alnST->queryStartUI =
      alnST->offsetUI
    + (endQueryAlnCStr - bestQueryCStr);

  switch(lastBitElmUC)
  { // Switch; check which sequence I am on the off base
    case defMoveDiagnol:   // Both sequences
      ++bestQueryCStr;
      ++bestRefCStr;
      break;

    case defMoveLeft:   // Last move was a deletion
      ++bestRefCStr;
      break;

    case defMoveUp:   // Last base was insertion
      ++bestQueryCStr;
  } // Switch; check which sequence I am on the off base

  /*******************************************************\
  * Fun-03 Sec-06 Sub-02: 
  *  - Add soft masking for the starting bases
  \*******************************************************/

  switch(softMaskBl)
  { // Switch: Check if applying softmaksing
    case 0: break;

    case 1:
      tmpQueryCStr = queryCStr;
      tmpRefCStr = refCStr;
      endUCPtr = alnST->alnAryUC + alnST->numAlnBasesUI;
    
      // Apply softmasking to the start region
      while(1)
      { // While I have softmasking to do
        if(tmpQueryCStr != bestQueryCStr)
        { // IF have a query base to mask
          *endUCPtr |= defSoftQueryFlag;
          ++tmpQueryCStr;
    
          if(tmpRefCStr != bestRefCStr)
          { // IF have a query base to mask
            *endUCPtr |= defSoftRefFlag;
            ++tmpRefCStr;
          } // IF have a query base to mask
        } // IF have a query base to mask
    
        else if(tmpRefCStr != bestRefCStr)
        { // IF have a query base to mask
          *endUCPtr |= defSoftRefFlag;
          ++tmpRefCStr;
        } // IF have a query base to mask
    
        else break; // Done
    
        ++endUCPtr;
        ++alnST->numBasesUI;
      } // While I have softmasking to do

      break;
    // Case 1: Apply softmaksing
  } // Switch: Check if applying softmaksing
  
  /*******************************************************\
  * Fun-03 Sec-06 Sub-02: 
  *  - clean up indels at end
  \*******************************************************/

  // This only happens when the gapopen penalty is 0
  switch(softMaskBl)
  { // Switch: converting hanging indels to softmasks?
    case 0: break;

    case 1:
    // Case 1: Conveting hanging indels to softmasks
      startUCPtr = alnST->alnAryUC;

      // Remove any hanging indels at the end
      while(!(*startUCPtr & 4))
      { // While I have no bases
         switch(*startUCPtr)
         { // Switch: Check wich soft mask I apply
           case defDelFlag:
             *startUCPtr = defSoftRefFlag;
             break;

           case defInsFlag:
             *startUCPtr = defSoftQueryFlag;
             break;
         } // Switch: Check wich soft mask I apply

         if(*startUCPtr == 0) break; // At end of the array

         ++startUCPtr;
      } // While I have no bases

      break;
    // Case 1: Conveting hanging indels to softmasks
  } // Switch: converting hanging indels to softmasks?

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Invert the error type array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  // flag the end of the error array
  startUCPtr = alnST->alnAryUC;
  endUCPtr = (alnST->alnAryUC + alnST->numBasesUI - 1);

  while(startUCPtr < endUCPtr)
  { // While I have elements to inver
    swapUC = *endUCPtr;
    *endUCPtr = *startUCPtr;
    *startUCPtr = swapUC;
    ++startUCPtr;
    --endUCPtr;
  } // While I have elements to inver

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-05:
  ^  - Add softmasking to the end and return
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  switch(softMaskBl)
  { // Switch: Check if doing softmasking for the end
    case 0: return alnST;

    case 1:
      // Move to the first base to sofmask
      ++endQueryAlnCStr;
      ++endRefAlnCStr;
      endUCPtr = alnST->alnAryUC + alnST->numBasesUI;

      // Apply softmasking to the start region
      while(1)
      { // While I have softmasking to do
        *endUCPtr = 0;

        if(*endQueryAlnCStr != '\0')
        { // IF have a query base to mask
          *endUCPtr |= defSoftQueryFlag;
          ++endQueryAlnCStr;

          if(*endRefAlnCStr != '\0')
          { // IF have a query base to mask
            *endUCPtr |= defSoftRefFlag;
            ++endRefAlnCStr;
            ++tmpRefCStr;
          } // IF have a query base to mask

          ++endUCPtr;
        } // IF have a query base to mask

        else if(*endRefAlnCStr != '\0')
        { // IF have a query base to mask
            *endUCPtr |= defSoftRefFlag;
            ++endRefAlnCStr;
            ++endUCPtr;
        } // IF have a query base to mask

        else  return alnST; // Done with softmaksing
      } // While I have softmasking to do

    default: return alnST;
  } // Switch: Check if doing softmasking for the end
} // cnvtDirMatrixToAlnAry

/*--------------------------------------------------------\
| Output:
|  - Prints
|    o Prints out the input alingment
\--------------------------------------------------------*/
void printAln(
  FILE *outFILE,      // File to print alingment to
  char *queryIdCStr,  // Id/name of query sequence
  char *queryAlnCStr, // Alinged query sequence
  char *refAlnCStr,   // Aligned reference sequence
  char *refIdCStr,    // Id/name of reference sequence
  long scoreL,        // Score of the alignment
  unsigned short lineWrapUS,
    // Number of characters per line (minimum is 42)
  struct alnStruct *alnST // alignment array in human form
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: printAln
   '  - Prints out the alignment
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Print out the header  for the alignment
   '  o fun-04 sec-04:
   '    - Print out tail of the alingment (missed by loop)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint8_t *alnAryUCPtr = 0;
   unsigned short wrapUS = lineWrapUS - 9;
     // - 9 to account for 9 spaces needed for seq entry
     // Min is 42 so that the comment section fits on one
     // one line
   uint32_t uiBase = 0; // Base priting out

   char headerCStr = "\
     \n# Eqx = Error line\
     \n#   - = is match\
     \n#   - X is mismatch\
     \n#   - I is insertion\
     \n#   - D is deletion\
     \n#   - S is soft mask on query and reference\
     \n#   - s is soft mask on query only\
     \n#   - P is soft mask on reference only\
     \n###########################################";

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Print out the header  for the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(wrapUS < 42) wrapUS = 42;

   fputs(
     "###########################################",
     outFILE
   );

   fprintf(outFILE, "# Query = %s", queryIdCStr);

   fprintf(
     outFILE,
     "#   - Query bases %u to %u",
     alnST->queryStartUI,
     alnST->queryEndUI
   );

   fprintf(outFILE, "# Ref = %s", refIdCStr);

   fprintf(
     outFILE,
     "#   - Reference bases %u to %u",
     alnST->refStartUI,
     alnST->refEndUI
   );
 
   fprintf(outFILE, "# Alignment Score = %l", scoreL);
   fputs(headerCStr, outFILE);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Print out the aligment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alnAryUCPtr = alnST->alnAryUC;
   uiBase += wrapUS; // Avoid overprinting with fwrite

   while(uiBase < alnST->numBasesUI)
   { // while I need to print the alignment
     fwrite("\n", sizeof(char), 1, outFILE);

     fwrite("Ref:     ", sizeof(char), 9, outFILE);
     fwrite(refAlnCStr, sizeof(char), wrapUS, outFILE);
     fwrite("\n", sizeof(char), 1, outFILE);

     fwrite("Query:   ", sizeof(char), 9, outFILE);
     fwrite(queryAlnCStr, sizeof(char), wrapUS, outFILE);
     fwrite("\n", sizeof(char), 1, outFILE);

     
     fwrite("Eqx:     ", sizeof(char), 9, outFILE);
     fwrite(alnAryUCPtr, sizeof(char), wrapUS, outFILE);
     fwrite("\n", sizeof(char), 1, outFILE);

     refAlnCStr += wrapUS;
     queryAlnCStr += wrapUS;
     alnAryUCPtr += wrapUS;
     uiBase += wrapUS;
   } // while I need to print the alignment

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-04:
   ^  - Print out tail of the alingment (missed by loop)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uiBase += 1 + wrapUS;

   if(uiBase < alnST->numBasesUI)
     if(uiBase > 1)
     { // If I have printed other parts of this alignment
       refAlnCStr += 1 - warpUS
       queryAlnCStr += 1 - warpUS
       alnAryUCPtr += 1 - warpUS
     } // If I have printed other parts of this alignment
     // Else; the alignment is covered in one wrap

     fwrite("\n", sizeof(char), 1, outFILE);

     fwrite("Ref:     ", sizeof(char), 9, outFILE);
     fwrite(refAlnCStr, sizeof(char), wrapUS, outFILE);
     fwrite("\n", sizeof(char), 1, outFILE);

     fwrite("Query:   ", sizeof(char), 9, outFILE);
     fwrite(queryAlnCStr, sizeof(char), wrapUS, outFILE);
     fwrite("\n", sizeof(char), 1, outFILE);

     
     fwrite("Eqx:     ", sizeof(char), 9, outFILE);
     fwrite(alnAryUCPtr, sizeof(char), wrapUS, outFILE);
     fwrite("\n", sizeof(char), 1, outFILE);
   } // If missed the last base

   return;
} // printAln

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: Sec-01: initAlnST
   '  - Initalize all values in alnST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   alnST->alnAryUC = 0;
   alnST->lenAlnAryUI = 0;
   alnST->numBasesUI = 0;
   alnST->refStartUI = 0;
   alnST->refEndUI = 0;
   alnST->queryStartUI = 0;
   alnST->queryEndUI = 0;

   return;
} // initAlnST

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-01: freeAlnST
   '  - Frees alnST and all variables in alnST
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(alnST->alnAryUC != 0) free(alnST->alnAryUC);

   if(heapBl != 0) free(alnST);

   return;
} // initAlnST
