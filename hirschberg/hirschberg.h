/*########################################################
# Name: hirschberg
# Use:
#  - Holds functions for doing a hirschberg global
#    alignment
# Libraries:
#  - "genScoreHirsch.h"               (No .c File)
#  o "genHirsch.h"                    (No .c File)
#  o "../general/alnStruct.h"         (No .c File)
#  o "../general/seqStruct.h"         (No .c File)
#  o "../general/alnSetStruct.h"      (No .c File)
#  o "../general/alnSeqDefaults.h"    (No .c File)
#  o "../general/genAln.h"            (No .c File)
#  o "../general/strToNum.h"          (No .c File)
#  o "../general/dataTypeShortHand.h" (No .c File)
#  o "../general/alnMatrixStruct.h"   (No .c File)
#  o "../general/genMath.h"           (No .c File)
# C Standard Libraries:
#  o <time.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <string.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - Includes and definitions
'  o fun-01 HirschbergFun:
'    - Does the recursive part of a Hirschberg alignment
'  o fun-02 Hirschberg:
'    - Sets up for and calls the recursvie function to
'      run a Hirschberg alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Includes and definitions
\-------------------------------------------------------*/

#ifndef HIRSCHBERG_H
#define HIRSCHBERG_H

#include "genScoreHirsch.h"

/*-------------------------------------------------------\
| Fun-01: HirschbergFun
|  - Do an Hirschberg alignment
| Input
|  - refSeqStr:
|    o reference sequence to align
|  - refStartUL:
|    o First base to align in the reference (index 0)
|  - refLenUL:
|    o Number of bases to align in the reference (index 1)
|  - qrySeqStr:
|    o query sequence to align
|  - qryStartUL:
|    o First base to align in the query (index 0)
|  - qryLenUL:
|    o Number of bases to align in the query (index 1)
|  - forScoreRowL:
|    o Row holding the forward (1st half of query) scores
|    o Must be the size of the full length reference
|  - revScoreRowL:
|    o Row holding the reverse (last half of query) scores
|    o Must be the size of the full length reference
|  - refAlnST:
|    o Holds the reference alignment and is a temporary
|      row for finding directions
|    o if -DHIRSCHTWOBIT is two bit array, else char array
|  - qryAlnST:
|    o Holds the query alignment
|    o if -DHIRSCHTWOBIT is two bit array, else char array
|  - dirRow:
|    o Is used for finding the directions. This is here to
|      ensure the Hirschberg is thread safe
|    o if -DHIRSCHTWOBIT is two bit array, else char array
|  - settings:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o refAlnST to hold the reference alignment
|    o qryAlnST to hold the query alignment
\-------------------------------------------------------*/
static void HirschbergFun(
  char *refSeqStr,  /*Reference sequence*/
  ulong refStartUL, /*1st reference base to align*/
  ulong refLenUL,   /*number of reference bases to align*/

  char *qrySeqStr,  /*Query sequence*/
  ulong qryStartUL, /*1st query base to align (index 0)*/
  ulong qryLenUL,   /*number of query bases to align*/

  long *forScoreRowL, /*Holds final forward row*/
  long *revScoreRowL, /*For finding reverse scores*/
  
  char *refAlnST,  /*Holds output reference alignment*/
  char *qryAlnST,  /*Holds the output query alignment*/
  char *dirRow,    /*Keeps thread safe*/

  struct alnSet *settings /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: HirschbergFun
   '  - Does the recursive part of a Hirschberg alignment
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Check if on a leaf (final part of alignment
   '  o fun-01 sec-03:
   '    - Get scores
   '  o fun-01 sec-04:
   '    - Find the midpoint
   '  o fun-01 sec-05:
   '    - Run the next hirschberg alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long forGapColL = 0;
   long revGapColL = 0;
   ulong midPointUL = 0;
   ulong ulFor = 0; /*Loop iterator (forward score)*/

   ulong ulGapOn = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Check if on a leaf (final part of alignment
   ^  o fun-01 sec-02 sub-01:
   ^    - Handle cases were I have just insertions
   ^  o fun-01 sec-02 sub-02:
   ^    - Handle cases were I have just deletions
   ^  o fun-01 sec-02 sub-03:
   ^    - Handle cases were I have to align last ref base
   ^  o fun-01 sec-02 sub-04:
   ^  - Handle cases were I have to align last query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Handle cases were I have just insertions
   \*****************************************************/

   if(refLenUL == 0)
   { /*If all remaing bases are query insertions*/
    for(
       ulGapOn = qryStartUL;
       ulGapOn < qryLenUL + qryStartUL;
       ++ulGapOn
     ) qryAlnST[ulGapOn] = defGapFlag;

     return; /*Nothing else to do*/
   } /*If all remaing bases are query insertions*/

   /*****************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Handle cases were I have just deletions
   \*****************************************************/

   if(qryLenUL == 0)
   { /*If all remaing bases are query deletions*/
    for(
       ulGapOn = refStartUL;
       ulGapOn < refLenUL + refStartUL;
       ++ulGapOn
    ) refAlnST[ulGapOn] = defGapFlag;

     return; /*Nothing else to do*/
   } /*If all remaing bases are query deletions*/

   /*****************************************************\
   * Fun-01 Sec-02 Sub-03:
   *  - Handle cases were I have to align last ref base
   \*****************************************************/

   if(refLenUL == 1)
   { /*If I have to align the last reference base*/
     if(qryLenUL == 0)
     { /*If bases are aligned (one ref & one query)*/
        refAlnST[refStartUL] = defGapFlag;
        return; /*Finished*/
     } /*If bases are aligned (one ref & one query)*/

     if(qryLenUL == 1)
     { /*If bases are aligned (one ref & one query)*/
        qryAlnST[qryStartUL] = defSnpFlag;
        refAlnST[refStartUL] = defSnpFlag;
        return; /*Finished*/
     } /*If bases are aligned (one ref & one query)*/

     positionSingleBase(
       *(refSeqStr + refStartUL),/*ref base*/
       refStartUL,             /*Position of ref base*/
       qrySeqStr,              /*first base of query*/
       qryStartUL,             /*positoin of query*/
       qryLenUL,               /*Length of the query*/
       refAlnST,               /*Array to hold alignment*/
       qryAlnST,               /*Array to hold alignment*/
       settings                /*Has Scoring variables*/
     );

     return; /*This base is now aligned*/
   } /*If I have to align the last reference base*/

   /*****************************************************\
   * Fun-01 Sec-02 Sub-04:
   *  - Handle cases were I have to align last query base
   \*****************************************************/

   if(qryLenUL == 1)
   { /*If I have to align the last query base*/
     if(refLenUL == 0)
     { /*If bases are aligned (one ref & one query)*/
        qryAlnST[qryStartUL] = defGapFlag;
        return; /*Finished*/
     } /*If bases are aligned (one ref & one query)*/

     if(refLenUL == 1)
     { /*If bases are aligned (one ref & one query)*/
        qryAlnST[qryStartUL] = defSnpFlag;
        refAlnST[refStartUL] = defSnpFlag;
        return; /*Finished*/
     } /*If bases are aligned (one ref & one query)*/

     positionSingleBase(
       *(qrySeqStr + qryStartUL),/*ref base*/
       qryStartUL,             /*Position of ref base*/
       refSeqStr,              /*first base of reference*/
       refStartUL,             /*positoin of query*/
       refLenUL,               /*Length of the query*/
       qryAlnST,               /*Array to hold alignment*/
       refAlnST,               /*Array to hold alignment*/
       settings                /*Has Scoring variables*/
     );

     return; /*Finshed aligning this query base*/
   } /*If I have to align the last query base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Get scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    forGapColL = 
      scoreHirschFor(
        refSeqStr,       /*Entire reference sequence*/
        refStartUL,      /*Starting base of ref target*/
        refLenUL,        /*length of ref target region*/
        qrySeqStr,       /*Query seq with coordinates*/
        qryStartUL,      /*Starting base of query target*/
        qryLenUL / 2,    /*Length of query target region*/
        forScoreRowL,    /*Array of scores to fill*/
        refAlnST,        /*direction row for gap extend*/
        settings         /*setttings to use*/
    ); /*Get the scores for the forward direction*/
    /*For -DNOGAPOPEN, refAlnST is ignored*/

    revGapColL = 
      scoreHirschRev(
        refSeqStr,     /*Entire reference sequence*/
        refStartUL,    /*Starting base of ref target*/
        refLenUL,      /*length of ref target region*/
        qrySeqStr,     /*Query seq with coordinates*/
        qryStartUL + (qryLenUL / 2),/*new query start*/
        qryLenUL - (qryLenUL / 2),  /*New query length*/
        revScoreRowL,  /*Array of scores to fill*/
        dirRow,        /*direction row for gap extend*/
        settings       /* setttings to use*/
      ); /*Get the scores for the reverse direction*/
      /*For -DNOGAPOPEN, dirRow is ignored*/
      /* I can get away with queryLen/2 here, because 
      `  queryLen is index 1 and the function takes in
      `  an lenth 1 argument
      `  I made this section thread safe by using refAlnST
      `    for the forward score and dirRow for the
      `    reverse row. So, you could run the forward and
      `    reverse scoring step in parrallel
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-04:
   ^   - Find the midpoint
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   forScoreRowL[refStartUL + refLenUL - 1] += revGapColL;
   midPointUL = refStartUL + refLenUL - 1;

   for(
      ulFor = refStartUL;
      ulFor < refStartUL + refLenUL - 1;
      ++ulFor
   ){ /*Loop; add up all scores*/
     forScoreRowL[ulFor] += revScoreRowL[ulFor + 1];
       /*The reverse row is already reversed*/

     if(forScoreRowL[ulFor] > forScoreRowL[midPointUL])
        midPointUL = ulFor;
   } /*Loop; add up all scores*/

   forGapColL += revScoreRowL[refStartUL];

   if(forGapColL > forScoreRowL[midPointUL])
      midPointUL = 0;
   else midPointUL = midPointUL + 1 - refStartUL;


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-05:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   HirschbergFun(
     refSeqStr,      /*Full reference sequence*/
     refStartUL,     /*Full reference sequence*/
     midPointUL,     /*Length of new reference sequence*/
     qrySeqStr,      /*Full query sequence*/
     qryStartUL,     /*Start of queyr target region*/
     qryLenUL / 2,   /*Length of query target region*/
     forScoreRowL,   /*For scoring*/
     revScoreRowL,   /*Has last line of scores*/
     refAlnST,       /*direction row for gap extend*/
     qryAlnST,       /*Holds the alignment codes*/
     dirRow,         /*For threadsafe scoreing*/
     settings        /*Settings for the alignment*/
   );

   HirschbergFun(
     refSeqStr,              /*Full reference sequence*/
     refStartUL + midPointUL,/*New reference start*/
     refLenUL - midPointUL,  /*New reference end*/
     qrySeqStr,              /*Full query sequence*/
     qryStartUL + (qryLenUL / 2),/*New query start*/
     qryLenUL - (qryLenUL / 2),  /*New query length*/
     forScoreRowL,           /*For scoring*/
     revScoreRowL,           /*Has last line of scores*/
     refAlnST,               /*Holds reference alingment*/
     qryAlnST,               /*Holds query alingment*/
     dirRow,                 /*For threadsafe scoring*/
     settings                /*Settings for alignment*/
   );

   return;
} /*HirschbergFun*/

/*-------------------------------------------------------\
| Fun-02: Hirschberg
|  - Sets up for and calls the recursvie function to
|    run a Hirschberg alignment
| Input:
|  - refST:
|    o Pointer to seqStruct structure with the reference
|      sequence and the first base (offsetUL) and last
|      base (endAlnUL) to align (both are index 0)
|  - qryST:
|    o Pointer to seqStruct structure with the query
|      sequence and the first base (offsetUL) and last
|      base (endAlnUL) to align (both are index 0)
|  - settings:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Returns:
|    o A alignment structure with the alignment
|    o 0 For memory errors
\-------------------------------------------------------*/
static struct alnStruct * Hirschberg(
  struct seqStruct *refST, /*Reference sequence to align*/
  struct seqStruct *qryST, /*Qeury sequence to align*/
  struct alnSet *settings  /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Hirschberg
   '  - Sets up for and calls the recursvie function to
   '    run a Hirschberg alignment
   '  o fun-02 sec-01:
   '    - Variable declerations
   '  o fun-02 sec-02:
   '    - Memory allocation (set up for Hirschberg)
   '  o fun-02 sec-03:
   '    - Run the hirschberg alignment
   '  o fun-02 sec-04:
   '   - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   ulong lenRefUL = refST->endAlnUL - refST->offsetUL + 1;
   ulong lenQryUL = qryST->endAlnUL - qryST->offsetUL + 1;
     /*+ 1 to convert to index 1 (values are index 0)*/

   long *forwardScoreRowL = 0;
   long *reverseScoreRowL = 0;

   struct alnStruct *alnST = 0;

   char *refAln = 0;
   char *qryAln = 0;
   char *dirRow = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-02:
   ^   - Memory allocation (set up for Hirschberg)
   ^   o fun-02 sec-02 sub-01:
   ^     - Initalize the ouput alignment structure 
   ^   o fun-02 sec-02 sub-02:
   ^     - Initalize the scoring rows
   ^   o fun-02 sec-02 sub-03:
   ^     - Initalize the direction rows
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-02 Sec-02 Sub-01:
   *  - Initalize the ouput alignment structure 
   \*****************************************************/

   refAln = calloc(lenRefUL + 1, sizeof(char));
   if(refAln == 0) return 0;
   *(refAln + lenRefUL) = defEndAlnFlag;

   dirRow = calloc(lenRefUL + 1, sizeof(char));

   if(dirRow == 0)
   { /*If I could not make another direction row*/
      free(refAln);
      refAln = 0;
      return 0;
   } /*If I could not make another direction row*/

   qryAln = calloc(lenQryUL + 1, sizeof(char));

   if(qryAln == 0)
   { /*If I could not make another direction row*/
      free(refAln);
      refAln = 0;

      free(dirRow);
      dirRow = 0;

      return 0;
   } /*If I could not make another direction row*/

   /*****************************************************\
   * Fun-02 Sec-02 Sub-02:
   *  - Initalize the scoring rows
   \*****************************************************/

   forwardScoreRowL = malloc(sizeof(long) * lenRefUL);

   if(forwardScoreRowL == 0)
   { /*If had a memory allocatoin error*/
     free(refAln);
     refAln = 0;

     free(qryAln);
     qryAln = 0;

     free(dirRow);
     dirRow = 0;

     return 0;
   } /*If had a memory allocatoin error*/

   reverseScoreRowL = malloc(sizeof(long) * lenRefUL);

   if(reverseScoreRowL == 0)
   { /*If had a memory allocatoin error*/
     free(refAln);
     refAln = 0;

     free(qryAln);
     qryAln = 0;

     free(dirRow);
     dirRow = 0;

     free(forwardScoreRowL);
     return 0;
   } /*If had a memory allocatoin error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-03:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Sending in offset values, because alignment array is
   ` sized to the alignmnet region
   */
   HirschbergFun(
     refST->seqCStr + refST->offsetUL,
     0,                /*1st reference base to align*/
     lenRefUL,         /*Length of ref region to align*/
     qryST->seqCStr + qryST->offsetUL,
     0,                /*1st query base to align*/
     lenQryUL,         /*length of query target region*/
     forwardScoreRowL, /*For scoring*/
     reverseScoreRowL, /*For scoring*/
     refAln,     /*Holds the reference alignment*/
     qryAln,     /*Holds the query alignment*/
     dirRow,     /*Direction row for thread safe scoring*/
     settings    /*Settings for the alignment*/
   );
     /*dirRow becomes a dummy variable for -DNOGAPOPEN*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-04:
   ^    - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(forwardScoreRowL);
   free(reverseScoreRowL);
 
   free(dirRow);
   dirRow = 0;


   alnST =
     hirschToAlnST(refST,qryST,settings,refAln,qryAln);

   free(refAln);
   refAln = 0;

   free(qryAln);
   qryAln = 0;

   return alnST; /*Is 0 if twoBitAlnToAlnSt failed*/
} /*Hirschberg*/

#endif
