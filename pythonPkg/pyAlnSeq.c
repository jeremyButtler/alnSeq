/*TODO: Figure out how to use PyArgs_ParseTupleAndKeywords
` TODO: Add in and improve mem-water
*/

/*#########################################################
# Name: pyAlnSeq
# Use:
#  - Holds wrapper functions to allow alnSeq to be used in
#    python
# Libraries:
#  o Most of the .h files in ../
# C Standard Libraries:
#  o <Python.h>
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>  // by alnSetStructure.h
#  o <string.h>
#  o <time.h>
#########################################################*/

/*Most of the template for this is from
` https://realpython.com/documenting-python-code/#documenting-your-python-code-base-using-docstrings
*/

#define PY_SSIZE_T_CLEAN /*s# in PyArg_ParseTuple*()*/
#include <Python.h>

#include "../hirschberg/hirschberg.h"
#include "../hirschberg/hirschbergNoGap.h"

#include "../memWater/memWater.h"
#include "../memWater/memWaterNoGap.h"

#include "../needleman/needleman.h"
#include "../needleman/needleNoGap.h"
#include "../needleman/needleTwoBit.h"
#include "../needleman/needleTwoBitNoGap.h"

#include "../waterman/waterman.h"
#include "../waterman/watermanNoGap.h"
#include "../waterman/waterTwoBit.h"
#include "../waterman/waterTwoBitNoGap.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' plAlnSeq SOF: Start Of Functions
'  - Python wrapper functions for alnSeq's pairwise
'    aligments
'  o fun-01 alnSeqHirsch:
'    - Runs an Hirschberg on two input sequences
'  o fun-02 alnSeqMemWater:
'    - Runs a memory effiecent Waterman alignment on the
'      input sequences (returns coordinates only)
'  o fun-03 alnSeqNeedle:
'    - Runs an Needleman Wunsch alignment on two input
'      sequences
'  o fun-04 alnSeqWater:
'    - Run an Waterman Smith alignment on two input
'      sequences
'  o struct-01 alnSeqFunST:
'    - Structer to hold functions names for linking
'  o struct-02 alnSeqModule:
'    - Structer to hold all function information for alnSeq
'  o init-01 PyInit_alnSeq:
'    - initializes the alnSeq extension for python
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: alnSeqHirsch (Fun-01:)
| Use:
|  - Runs an Hirschberg on two input sequences
| Input:
|  - ref:
|    o Reference sequence
|  - query:
|    o query sequence
|  - gapopen:
|    o Gap opening penalty
|  - gapextend:
|    o Gap extension penalty
|  - refStart:
|    o Starting position to align reference at
|  - refEnd:
|    o Position to stop aligning reference at
|  - queryStart:
|    o Starting position to align queryerence at
|  - queryEnd:
|    o Position to stop aligning query at
|  - scoreMatrix 
|    o Path to scoring matrix for a the alignment
| Output:
|  - Returns:
|    o python tuple with the aligned reference sequence
|      (1st position) and query sequence (2nd position)
\--------------------------------------------------------*/
static PyObject * alnSeqHirsch(
   PyObject *self,
   PyObject *args, /*Arguments from user*/
   PyObject *kw    /*Key words to get agruments with*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC:
   '  - Run an Hirschberg alignment on two sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Get and check user input
   '  o fun-01 sec-03:
   '    - Get read lengths and set up for alignment
   '  o fun-01 sec-04:
   '    - Run the Hirschberg
   '  o fun-01 sec-05:
   '    - Convert the alignment to strings
   '  o fun-01 sec-06:
   '    - Convert the aligned c-strings to python strings
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *scoreMatrixStr = 0;
   char *matchMatrixStr = 0;
   char *alnRefStr = 0;
   char *alnQryStr = 0;
   long errUL = 0;

   struct seqStruct refST;
   struct seqStruct qryST;

   struct alnSet settings;
   struct alnStruct *alnST = 0;

   FILE *scoreFILE = 0;
   FILE *matchFILE = 0;

   PyObject *retObj;

   char * keywordsAry[] =
      {
         "ref",
         "query",
         "gapOpen",
         "gapExtend",
         "refStart",
         "refEnd",
         "queryStart",
         "queryEnd",
         "noGapBool",
         "scoreMatrix",
         "matchMatrix",
         NULL
      }; /*Keywords to look up user input with*/
         
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Get and check user input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&settings);
   initSeqST(&refST);
   initSeqST(&qryST);

   /*Get the user input*/
   errUL =
      PyArg_ParseTupleAndKeywords(
         args,
         kw,
         "s#s#|hhkkkkBss",
         keywordsAry,

         &(refST.seqCStr),
         &(refST.lenSeqUL),
         &(qryST.seqCStr),
         &(qryST.lenSeqUL),

         &settings.gapOpenC,
         &settings.gapExtendC,

         &(refST.offsetUL),
         &(refST.endAlnUL),
         &(qryST.offsetUL),
         &(qryST.endAlnUL),

         &(settings.noGapBl),

         &scoreMatrixStr,
         &matchMatrixStr
      );

   if(!errUL)
   { /*If: no user input was input*/
      if(refST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: No Reference sequence input\n"
         );
      else if(qryST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: No query sequence input\n"
         );
      else
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Invalid Input\n"
         );

      return NULL;
   } /*If: no user input was input*/

   if(scoreMatrixStr != 0)
   { /*If: the user provided a matrix*/
      scoreFILE = fopen(scoreMatrixStr, "r");

      if(scoreFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Unable to open scoring matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInScoreFile(&settings, scoreFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Scoring file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   if(matchMatrixStr != 0)
   { /*If: the user provided a matrix*/
      matchFILE = fopen(matchMatrixStr, "r");

      if(matchFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Unable to open match matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInMatchFile(&settings, matchFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Match file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Get read lengths and set up for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(refST.endAlnUL == 0)
      refST.endAlnUL = refST.lenSeqUL - 1;
   if(qryST.endAlnUL == 0)
      qryST.endAlnUL = qryST.lenSeqUL - 1;

   if(settings.noGapBl != 0) settings.noGapBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Run the hirschberg
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*This speeds things up*/
   seqToLookupIndex(refST.seqCStr);
   seqToLookupIndex(qryST.seqCStr);

   if(settings.noGapBl)
      alnST = HirschbergNoGap(&refST, &qryST, &settings);
   else
      alnST = Hirschberg(&refST, &qryST, &settings);

   lookupIndexToSeq(qryST.seqCStr);
   lookupIndexToSeq(refST.seqCStr);

   if(alnST == 0)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-05:
   ^  - Convert the alignment to strings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errUL =
     alnSTToSeq(
        &refST, /*Has reference sequence*/
        &qryST, /*Has the query sequence*/
        alnST,  /*Has the alignment*/
        &settings, /*Keep the entire alignment*/
        &alnRefStr,/*Will hold reference aligned sequence*/
        &alnQryStr /*Will hold the query aligned sequence*/
   ); /*Convert an alnStruct to two aligned sequences*/

   freeAlnST(alnST);
   alnST = 0;

   if(errUL)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-06:
   ^  - Convert the aligned c-strings to python strings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Convert aligned sequence c-strs to python strings*/
   /*This converts it to a tuple*/
   retObj =
     Py_BuildValue("[ssl]",alnRefStr,alnQryStr,0);

   /*These frees may mess up my code. Not sure*/
   free(alnRefStr);
   alnRefStr = 0;

   free(alnQryStr);
   alnQryStr = 0;

   return retObj;
} /*alnSeqHirsch*/

/*--------------------------------------------------------\
| Name: alnSeqMemWater (Fun-02:)
| Use:
|  - Runs an Memory efficent Waterman on two sequences
| Input:
|  - ref:
|    o Reference sequence
|  - query:
|    o query sequence
|  - gapopen:
|    o Gap opening penalty
|  - gapextend:
|    o Gap extension penalty
|  - refStart:
|    o Starting position to align reference at
|  - refEnd:
|    o Position to stop aligning reference at
|  - queryStart:
|    o Starting position to align queryerence at
|  - queryEnd:
|    o Position to stop aligning query at
|  - scoreMatrix 
|    o Path to scoring matrix for a the alignment
|  - noGapBool:
|    o true: Only use gap opening penalties; no extensions
|    o false: Use gap extension penalties
| Output:
|  - Returns:
|    o python tuple with the first reference base (1st)
|      last reference base (2nd), first query base (3rd),
|      last query base (4th), and score (5th)
\--------------------------------------------------------*/
static PyObject * alnSeqMemWater(
   PyObject *self,
   PyObject *args, /*Arguments from user*/
   PyObject *kw    /*Key words to get agruments with*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC:
   '  - Run an memory efficent Waterman on two sequences
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Get and check user input
   '  o fun-01 sec-03:
   '    - Get read lengths and set up for alignment
   '  o fun-01 sec-04:
   '    - Run the memWater aligner
   '  o fun-01 sec-05:
   '    - Convert the alignment to strings
   '  o fun-01 sec-06:
   '    - Convert the aligned c-strings to python strings
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *scoreMatrixStr = 0;
   char *matchMatrixStr = 0;
   long errUL = 0;

   struct seqStruct refST;
   struct seqStruct qryST;

   struct alnSet settings;
   struct alnMatrix *matrixST;

   ulong refStartUL = 0;
   ulong refEndUL = 0;
   ulong qryStartUL = 0;
   ulong qryEndUL = 0;

   FILE *scoreFILE = 0;
   FILE *matchFILE = 0;

   PyObject *retObj;

   char * keywordsAry[] =
      {
         "ref",
         "query",
         "gapOpen",
         "gapExtend",
         "refStart",
         "refEnd",
         "queryStart",
         "queryEnd",
         "noGapBool",
         "scoreMatrix",
         "matchMatrix",
         NULL
      }; /*Keywords to look up user input with*/
         
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-02:
   ^  - Get and check user input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&settings);
   initSeqST(&refST);
   initSeqST(&qryST);

   /*Get the user input*/
   errUL =
      PyArg_ParseTupleAndKeywords(
         args,
         kw,
         "s#s#|hhkkkkBss",
         keywordsAry,
         &(refST.seqCStr),
         &(refST.lenSeqUL),
         &(qryST.seqCStr),
         &(qryST.lenSeqUL),

         &settings.gapOpenC,
         &settings.gapExtendC,

         &(refST.offsetUL),
         &(refST.endAlnUL),
         &(qryST.offsetUL),
         &(qryST.endAlnUL),

         &(settings.noGapBl),

         &scoreMatrixStr,
         &matchMatrixStr
      );

   if(!errUL)
   { /*If: no user input was input*/
      if(refST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: No Reference sequence input\n"
         );
      else if(qryST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: No query sequence input\n"
         );
      else
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Invalid Input\n"
         );

      return NULL;
   } /*If: no user input was input*/

   if(scoreMatrixStr != 0)
   { /*If: the user provided a matrix*/
      scoreFILE = fopen(scoreMatrixStr, "r");

      if(scoreFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Unable to open scoring matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInScoreFile(&settings, scoreFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Scoring file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   if(matchMatrixStr != 0)
   { /*If: the user provided a matrix*/
      matchFILE = fopen(matchMatrixStr, "r");

      if(matchFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Unable to open match matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInMatchFile(&settings, matchFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Match file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-03:
   ^  - Get read lengths and set up for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(refST.endAlnUL == 0)
      refST.endAlnUL = refST.lenSeqUL - 1;
   if(qryST.endAlnUL == 0)
      qryST.endAlnUL = qryST.lenSeqUL - 1;

   if(settings.noGapBl != 0) settings.noGapBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-04:
   ^  - Run the memory efficent Waterman
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*This speeds things up*/
   seqToLookupIndex(refST.seqCStr);
   seqToLookupIndex(qryST.seqCStr);

   if(settings.noGapBl)
      matrixST = memWaterNoGap(&qryST, &refST, &settings);
   else
      matrixST = memWater(&qryST, &refST, &settings);

   lookupIndexToSeq(refST.seqCStr);
   lookupIndexToSeq(qryST.seqCStr);

   if(matrixST == 0)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-01 Sec-06:
   ^  - Convert the aligned c-strings to python strings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   indexToCoord(
      matrixST->lenRefUL,
      matrixST->bestStartIndexUL,
      refStartUL,
      qryStartUL
   );

   indexToCoord(
      matrixST->lenRefUL,
      matrixST->bestEndIndexUL,
      refEndUL,
      qryEndUL
   );

   /*Convert aligned sequence c-strs to python strings*/
   /*This converts it to a tuple*/
   retObj =
     Py_BuildValue(
        "[lkkkk]",
        matrixST->bestScoreL,
        refStartUL,
        refEndUL,
        qryStartUL,
        qryEndUL
   );

   /*These frees may mess up my code. Not sure*/
   freeAlnMatrix(matrixST);
   matrixST = 0;

   return retObj;
} /*alnSeqMemWater*/

/*--------------------------------------------------------\
| Name: alnSeqNeedle (Fun-03:)
| Use:
|  - Runs an Needleman Wunsch alignment on two input
|    sequences
| Input:
|  - ref:
|    o Reference sequence
|  - query:
|    o query sequence
|  - gapopen:
|    o Gap opening penalty
|  - gapextend:
|    o Gap extension penalty
|  - refStart:
|    o Starting position to align reference at
|  - refEnd:
|    o Position to stop aligning reference at
|  - queryStart:
|    o Starting position to align queryerence at
|  - queryEnd:
|    o Position to stop aligning query at
|  - scoreMatrix 
|    o Path to scoring matrix for a the alignment
| Output:
|  - Returns:
|    o python tuple with the aligned reference sequence
|      (1st position), query sequence (2nd position), and
|      the score (3rd position)
\--------------------------------------------------------*/
static PyObject * alnSeqNeedle(
   PyObject *self,
   PyObject *args, /*Arguments from user*/
   PyObject *kw    /*Key words to get agruments with*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC:
   '  - Run a Needleman alignment on two sequences
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Get and check user input
   '  o fun-03 sec-03:
   '    - Get read lengths and set up for alignment
   '  o fun-03 sec-04:
   '    - Run the Needleman
   '  o fun-03 sec-05:
   '    - Convert the alignment to strings
   '  o fun-03 sec-06:
   '    - Convert the aligned c-strings to python strings
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *scoreMatrixStr = 0;
   char *matchMatrixStr = 0;
   char *alnRefStr = 0;
   char *alnQryStr = 0;
   long errUL = 0;
   long bestScoreL = 0;

   struct seqStruct refST;
   struct seqStruct qryST;

   struct alnMatrix *alnMtrxST = 0;
   struct alnMatrixTwoBit *alnMtrxTBST = 0;
   struct alnSet settings;
   struct alnStruct *alnST = 0;

   FILE *scoreFILE = 0;
   FILE *matchFILE = 0;

   PyObject *retObj;

   char * keywordsAry[] =
      {
         "ref",
         "query",
         "gapOpen",
         "gapExtend",
         "refStart",
         "refEnd",
         "queryStart",
         "queryEnd",
         "noGapBool",
         "twoBitBool",
         "scoreMatrix",
         "matchMatrix",
         NULL
      }; /*Keywords to look up user input with*/
         
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Get and check user input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&settings);
   initSeqST(&refST);
   initSeqST(&qryST);

   /*Get the user input*/
   errUL =
      PyArg_ParseTupleAndKeywords(
         args,
         kw,
         "s#s#|hhkkkkBBss",
         keywordsAry,
         &(refST.seqCStr),
         &(refST.lenSeqUL),
         &(qryST.seqCStr),
         &(qryST.lenSeqUL),

         &settings.gapOpenC,
         &settings.gapExtendC,

         &(refST.offsetUL),
         &(refST.endAlnUL),
         &(qryST.offsetUL),
         &(qryST.endAlnUL),

         &(settings.noGapBl),
         &(settings.twoBitBl),

         &scoreMatrixStr,
         &matchMatrixStr
      );

   if(!errUL)
   { /*If: no user input was input*/
      if(refST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqeNeedle: No Reference sequence input\n"
         );
      else if(qryST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqNeedle: No query sequence input\n"
         );
      else
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Invalid Input\n"
         );

      return NULL;
   } /*If: no user input was input*/

   if(scoreMatrixStr != 0)
   { /*If: the user provided a matrix*/
      scoreFILE = fopen(scoreMatrixStr, "r");

      if(scoreFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqNeedle: Unable to open scoring matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInScoreFile(&settings, scoreFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqdNeedle: Scoring file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   if(matchMatrixStr != 0)
   { /*If: the user provided a matrix*/
      matchFILE = fopen(matchMatrixStr, "r");

      if(matchFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Unable to open match matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInMatchFile(&settings, matchFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Match file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Get read lengths and set up for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(refST.endAlnUL == 0)
      refST.endAlnUL = refST.lenSeqUL - 1;
   if(qryST.endAlnUL == 0)
      qryST.endAlnUL = qryST.lenSeqUL - 1;

   if(settings.noGapBl != 0) settings.noGapBl = 1;
   if(settings.twoBitBl != 0) settings.twoBitBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-04:
   ^  - Run the Needleman
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*This speeds things up*/
   seqToLookupIndex(refST.seqCStr);
   seqToLookupIndex(qryST.seqCStr);

   /*Get the directional matrix*/
   if(settings.noGapBl && settings.twoBitBl)
      alnMtrxTBST = 
         NeedleTwoBitNoGap(&qryST,&refST,&settings);
   else if(settings.twoBitBl)
      alnMtrxTBST = NeedleTwoBit(&qryST,&refST,&settings);
   else if(settings.noGapBl)
      alnMtrxST =NeedleAlnNoGap(&qryST, &refST, &settings);
   else
      alnMtrxST = NeedlemanAln(&qryST, &refST, &settings);

   lookupIndexToSeq(refST.seqCStr);
   lookupIndexToSeq(qryST.seqCStr);

   if(alnMtrxST == 0 && alnMtrxTBST == 0)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*convert directional matrix to an alignment struct*/

   if(settings.twoBitBl)
   { /*If: I used two bit arrays*/
      alnST =
         twoBitDirMatrixToAln(
            &refST,
            &qryST,
            alnMtrxTBST->bestEndIndexUL,
            &settings,
            alnMtrxTBST
      );

      bestScoreL = alnMtrxTBST->bestScoreL;
      freeAlnMatrixTwoBit(alnMtrxTBST);
      alnMtrxST = 0;
   } /*If: I used two bit arrays*/

   else
   { /*Else: I used byte arrays*/
      alnST =
         dirMatrixToAln(
            &refST,
            &qryST,
            alnMtrxST->bestEndIndexUL,
            &settings,
            alnMtrxST
      );

      bestScoreL = alnMtrxST->bestScoreL;
      freeAlnMatrix(alnMtrxST); /*No longer need*/
      alnMtrxST = 0;
   } /*Else: I used byte arrays*/

   if(alnST == 0)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-05:
   ^  - Convert the alignment to strings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errUL =
     alnSTToSeq(
        &refST, /*Has reference sequence*/
        &qryST, /*Has the query sequence*/
        alnST,  /*Has the alignment*/
        &settings, /*Keep the entire alignment*/
        &alnRefStr,/*Will hold reference aligned sequence*/
        &alnQryStr /*Will hold the query aligned sequence*/
   ); /*Convert an alnStruct to two aligned sequences*/

   freeAlnST(alnST);
   alnST = 0;

   if(errUL)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-06:
   ^  - Convert the aligned c-strings to python strings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Convert aligned sequence c-strs to python strings*/
   /*This converts it to a tuple*/
   retObj =
     Py_BuildValue("[ssl]",alnRefStr,alnQryStr,bestScoreL);

   /*These frees may mess up my code. Not sure*/
   free(alnRefStr);
   alnRefStr = 0;

   free(alnQryStr);
   alnQryStr = 0;

   return retObj;
} /*alnSeqNeedle*/

/*--------------------------------------------------------\
| Name: alnSeqWater (Fun-04:)
| Use:
|  - Run an Waterman Smith alignment on two input sequences
| Input:
|  - ref:
|    o Reference sequence
|  - query:
|    o query sequence
|  - gapopen:
|    o Gap opening penalty
|  - gapextend:
|    o Gap extension penalty
|  - refStart:
|    o Starting position to align reference at
|  - refEnd:
|    o Position to stop aligning reference at
|  - queryStart:
|    o Starting position to align queryerence at
|  - queryEnd:
|    o Position to stop aligning query at
|  - scoreMatrix 
|    o Path to scoring matrix for a the alignment
|  - fullAln:
|    o Print out the full alignment.
| Output:
|  - Returns:
|    o python tuple with the aligned reference sequence
|      (1st position), query sequence (2nd position), and
|      the score (3rd position)
\--------------------------------------------------------*/
static PyObject * alnSeqWater(
   PyObject *self,
   PyObject *args, /*Arguments from user*/
   PyObject *kw    /*Key words to get agruments with*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC:
   '  - Run an Waterman alignment on two sequences
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-01:
   '    - Get and check user input
   '  o fun-04 sec-03:
   '    - Get read lengths and set up for alignment
   '  o fun-04 sec-04:
   '    - Run the Waterman
   '  o fun-04 sec-05:
   '    - Convert the alignment to strings
   '  o fun-04 sec-06:
   '    - Convert the aligned c-strings to python strings
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *scoreMatrixStr = 0;
   char *matchMatrixStr = 0;
   char *alnRefStr = 0;
   char *alnQryStr = 0;
   long errUL = 0;
   long bestScoreL = 0;

   struct seqStruct refST;
   struct seqStruct qryST;

   struct alnMatrix *alnMtrxST = 0;
   struct alnMatrixTwoBit *alnMtrxTBST = 0;
   struct alnSet settings;
   struct alnStruct *alnST = 0;

   FILE *scoreFILE = 0;
   FILE *matchFILE = 0;

   PyObject *retObj;

   char * keywordsAry[] =
      {
         "ref",
         "query",
         "gapOpen",
         "gapExtend",
         "refStart",
         "refEnd",
         "queryStart",
         "queryEnd",
         "fullAln",
         "noGapBool",
         "twoBitBool",
         "scoreMatrix",
         "matchMatrix",
         NULL /*Needed as end*/
      }; /*Keywords to look up user input with*/
         
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Get and check user input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&settings);
   initSeqST(&refST);
   initSeqST(&qryST);

   /*Get the user input*/
   errUL =
      PyArg_ParseTupleAndKeywords(
         args,
         kw,
         "s#s#|hhkkkkBBBss",
         keywordsAry,

         &(refST.seqCStr),
         &(refST.lenSeqUL),
         &(qryST.seqCStr),
         &(qryST.lenSeqUL),

         &settings.gapOpenC,
         &settings.gapExtendC,

         &(refST.offsetUL),
         &(refST.endAlnUL),
         &(qryST.offsetUL),
         &(qryST.endAlnUL),

         &(settings.pFullAlnBl),
         &(settings.noGapBl),
         &(settings.twoBitBl),

         &scoreMatrixStr,
         &matchMatrixStr
      );

   if(!errUL)
   { /*If: no user input was input*/
      if(refST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqWater: No Reference sequence input\n"
         );
      else if(qryST.seqCStr == 0)
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqWater: No query sequence input\n"
         );
      else
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqWater: Invalid Input\n"
         );

      return NULL;
   } /*If: no user input was input*/

   if(scoreMatrixStr != 0)
   { /*If: the user provided a matrix*/
      scoreFILE = fopen(scoreMatrixStr, "r");

      if(scoreFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqWater: Unable to open scoring matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInScoreFile(&settings, scoreFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqWater: Scoring file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   if(matchMatrixStr != 0)
   { /*If: the user provided a matrix*/
      matchFILE = fopen(matchMatrixStr, "r");

      if(matchFILE == 0)
      { /*If: the scoring matrix file was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Unable to open match matrix\n"
         );
         return NULL;
      } /*If: the scoring matrix file was invalid*/

      errUL = readInMatchFile(&settings, matchFILE);

      if(errUL)
      { /*If: the scoring matrix was invalid*/
         PyErr_SetString(
            PyExc_ValueError,
            "alnSeqHirsch: Match file is invalid\n"
         );
         return NULL;
      } /*If: the scoring matrix was invalid*/
   } /*If: the user provided a matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Get read lengths and set up for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(refST.endAlnUL == 0)
      refST.endAlnUL = refST.lenSeqUL - 1;
   if(qryST.endAlnUL == 0)
      qryST.endAlnUL = qryST.lenSeqUL - 1;

   if(settings.pFullAlnBl != 0) settings.pFullAlnBl = 1;
   if(settings.noGapBl != 0) settings.noGapBl = 1;
   if(settings.twoBitBl != 0) settings.twoBitBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-04:
   ^  - Run the Waterman
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*This speeds things up*/
   seqToLookupIndex(refST.seqCStr);
   seqToLookupIndex(qryST.seqCStr);

   /*Get the directional matrix*/

   if(settings.noGapBl && settings.twoBitBl)
      alnMtrxTBST = 
         WaterTwoBitNoGap(&qryST,&refST,&settings);

   else if(settings.twoBitBl)
      alnMtrxTBST = WaterTwoBit(&qryST,&refST,&settings);

   else if(settings.noGapBl)
      alnMtrxST =
         WatermanAlnNoGap(&qryST, &refST, &settings);
   else
      alnMtrxST = WatermanAln(&qryST, &refST, &settings);

   lookupIndexToSeq(refST.seqCStr);
   lookupIndexToSeq(qryST.seqCStr);

   if(alnMtrxST == 0 && alnMtrxTBST == 0)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*convert directional matrix to an alignment struct*/

   if(settings.twoBitBl)
   { /*If: I used two bit arrays*/
      alnST =
         twoBitDirMatrixToAln(
            &refST,
            &qryST,
            alnMtrxTBST->bestEndIndexUL,
            &settings,
            alnMtrxTBST
      );

      bestScoreL = alnMtrxTBST->bestScoreL;
      freeAlnMatrixTwoBit(alnMtrxTBST);
      alnMtrxST = 0;
   } /*If: I used two bit arrays*/

   else
   { /*Else: I used byte arrays*/
      alnST =
         dirMatrixToAln(
            &refST,
            &qryST,
            alnMtrxST->bestEndIndexUL,
            &settings,
            alnMtrxST
      );

      bestScoreL = alnMtrxST->bestScoreL;
      freeAlnMatrix(alnMtrxST); /*No longer need*/
      alnMtrxST = 0;
   } /*Else: I used byte arrays*/

   if(alnST == 0)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-05:
   ^  - Convert the alignment to strings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errUL =
     alnSTToSeq(
        &refST, /*Has reference sequence*/
        &qryST, /*Has the query sequence*/
        alnST,  /*Has the alignment*/
        &settings,
        &alnRefStr,/*holds reference aligned sequence*/
        &alnQryStr /*holds the query aligned sequence*/
   ); /*Convert an alnStruct to two aligned sequences*/

   freeAlnST(alnST);
   alnST = 0;

   if(errUL)
   { /*If: Had a memory error*/
      PyErr_NoMemory();
      return NULL;
   } /*If: Had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-06:
   ^  - Convert the aligned c-strings to python strings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Convert aligned sequence c-strs to python strings*/
   /*This converts it to a tuple*/
   retObj =
     Py_BuildValue("[ssl]",alnRefStr,alnQryStr,bestScoreL);

   /*These frees may mess up my code. Not sure*/
   free(alnRefStr);
   alnRefStr = 0;

   free(alnQryStr);
   alnQryStr = 0;

   return retObj;
} /*alnSeqWater*/

/*--------------------------------------------------------\
| Name: alnSeqFunST [Struct-01:]
| Use:
|  - Structer to hold functions names for linking
|  o struct-01 sec-01:
|    - Set up for Hirschberg
|  o struct-01 sec-02:
|    - Set up for the memory efficent Smith Waterman
|  o struct-01 sec-03:
|    - Set up for the Needleman Wunsch
|  o struct-01 sec-04:
|    - Set up for the Waterman Smith
\--------------------------------------------------------*/
static PyMethodDef alnSeqFunST[] =
{

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Struct-01 Sec-01:
   ^  - Set up for Hirschberg
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  { /*Hirschberg wrapper*/
     "alnSeqHirsch",              /*Python function name*/
     alnSeqHirsch,                /*C function name*/
     METH_VARARGS | METH_KEYWORDS,/*takes keywords & args*/
     "Runs a Hirschberg on the input reference sequence\n\
        and query sequence.\n\
      This functions returns an list with the aligned\n\
        reference sequence (list[0]) and the aligned\n\
        query sequence (list[1])\n\
      alnSeqHirsch(ref=ref, query=query, [options...])\n\
        Options:\n\
          ref: [ref = sequence; Required]\n\
            - Reference sequence to align.\n\
          Query: [query = sequence; Required]\n\
            - Query sequence to align.\n\
          gapOpen: [gapOpen = -10]\n\
            - Score for starting an gap (indel).\n\
          gapExtend: [gapExtend = -1]\n\
            - Score for extending an gap (indel).\n\
          refStart: [refStart = 0]\n\
            - First base to align in reference.\n\
          refEnd: [refEnd = length(reference) -  1]\n\
            - Last base to align in reference (index 0)\n\
          queryStart: [queryStart = 0]\n\
            - First base to align in query.\n\
          queryEnd: [queryEnd = length(reference) - 1]\n\
            - Last base to align in query (index 0).\n\
          noGapBool: [false]\n\
            - Do not use gap extension penalties\n\
          scoreMatrix: [scoreMatrix = NULL]\n\
            - file name with scoring matrix to use.\n\
            - Default is the EDNAFULL matrix\n\
          matchMatrix: [matchMatrix = NULL]\n\
            - file name with match matrix to use.\n\
            - Default is DNA\n\
     " /*Documentation*/
  }, /*Hirschberg wrapper*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Struct-01 Sec-02:
   ^  - Set up for memWater
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  { /*memWater wrapper*/
     "alnSeqMemWater",            /*Python function name*/
     alnSeqMemWater,              /*C function name*/
     METH_VARARGS | METH_KEYWORDS,/*takes keywords & args*/
     "Runs a Memory efficent Waterman on the input\n\
      reference sequence and query sequence.\n\
      This functions returns the starting and ending\n\
        coordinates and the score as a list\n\
      alnSeqMemWater(ref=ref,query=query,[options...])\n\
        Options:\n\
          ref: [ref = sequence; Required]\n\
            - Reference sequence to align.\n\
          Query: [query = sequence; Required]\n\
            - Query sequence to align.\n\
          gapOpen: [gapOpen = -10]\n\
            - Score for starting an gap (indel).\n\
          gapExtend: [gapExtend = -1]\n\
            - Score for extending an gap (indel).\n\
          refStart: [refStart = 0]\n\
            - First base to align in reference.\n\
          refEnd: [refEnd = length(reference) -  1]\n\
            - Last base to align in reference (index 0)\n\
          queryStart: [queryStart = 0]\n\
            - First base to align in query.\n\
          queryEnd: [queryEnd = length(reference) - 1]\n\
            - Last base to align in query (index 0).\n\
          noGapBool: [false]\n\
            - Do not use gap extension penalties\n\
          scoreMatrix: [scoreMatrix = NULL]\n\
            - file name with scoring matrix to use.\n\
            - Default is the EDNAFULL matrix\n\
          matchMatrix: [matchMatrix = NULL]\n\
            - file name with match matrix to use.\n\
            - Default is DNA\n\
     " /*Documentation*/
  }, /*MemWater wrapper*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Struct-01 Sec-03:
   ^  - Set up for the Needleman Wunsch
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   {
     "alnSeqNeedle",
     alnSeqNeedle,
     METH_VARARGS | METH_KEYWORDS,/*takes keywords & args*/
     "Runs an Needleman on the input reference sequence\n\
        and query sequence.\n\
      This functions returns an list with the aligned\n\
        reference sequence (list[0]), the aligned\n\
        query sequence (list[1]), and the score for the\n\
        alignment (list[2])\n\
      alnSeqNeedle(ref=ref, query=query, [options...])\n\
        Options:\n\
          ref: [ref = sequence; Required]\n\
            - Reference sequence to align.\n\
          Query: [query = sequence; Required]\n\
            - Query sequence to align.\n\
          gapOpen: [gapOpen = -10]\n\
            - Score for starting an gap (indel).\n\
          gapExtend: [gapExtend = -1]\n\
            - Score for extending an gap (indel).\n\
          refStart: [refStart = 0]\n\
            - First base to align in reference.\n\
          refEnd: [refEnd = length(reference) -  1]\n\
            - Last base to align in reference (index 0)\n\
          queryStart: [queryStart = 0]\n\
            - First base to align in query.\n\
          queryEnd: [queryEnd = length(reference) - 1]\n\
            - Last base to align in query (index 0).\n\
          noGapBool: [false]\n\
            - Do not use gap extension penalties\n\
          twoBitBool: [false]\n\
            - Use two bit arrays, which are slower\n\
              (2x), but use less memory\n\
          scoreMatrix: [scoreMatrix = NULL]\n\
            - file name with scoring matrix to use.\n\
            - Default is the EDNAFULL matrix\n\
          matchMatrix: [matchMatrix = NULL]\n\
            - file name with match matrix to use.\n\
            - Default is DNA\n\
     "
   },

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Struct-01 Sec-04:
   ^  - Set up for the Waterman Smith
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   {
     "alnSeqWater",
     alnSeqWater,
     METH_VARARGS | METH_KEYWORDS,
     "Runs an Waterman alignment on the input reference\n\
        sequence and query sequence.\n\
      This functions returns an list with the aligned\n\
        reference sequence (list[0]), the aligned\n\
        query sequence (list[1]), and the score for the\n\
        alignment (list[2])\n\
      alnSeqWater(ref=ref, query=query, [options...])\n\
        Options:\n\
          ref: [ref = sequence; Required]\n\
            - Reference sequence to align.\n\
          Query: [query = sequence; Required]\n\
            - Query sequence to align.\n\
          gapOpen: [gapOpen = -10]\n\
            - Score for starting an gap (indel).\n\
          gapExtend: [gapExtend = -1]\n\
            - Score for extending an gap (indel).\n\
          refStart: [refStart = 0]\n\
            - First base to align in reference.\n\
          refEnd: [refEnd = length(reference) -  1]\n\
            - Last base to align in reference (index 0).\n\
          queryStart: [queryStart = 0]\n\
            - First base to align in query.\n\
          queryEnd: [queryEnd = length(reference) - 1]\n\
            - Last base to align in query (index 0).\n\
          fullAln: [fullAln = False]\n\
            - Print out the full alignment instead of\n\
              just printing out the aligned region.\n\
          noGapBool: [false]\n\
            - Do not use gap extension penalties\n\
          twoBitBool: [false]\n\
            - Use two bit arrays, which are slower\n\
              (2x), but use less memory\n\
          scoreMatrix: [scoreMatrix = NULL]\n\
            - file name with scoring matrix to use.\n\
            - Default is the EDNAFULL matrix\n\
          matchMatrix: [matchMatrix = NULL]\n\
            - file name with match matrix to use.\n\
            - Default is DNA\n\
     "
  },

  {NULL,NULL,0,NULL}
     /*Need a null item , otherwise python import will
     ` will error out on the initialize import function
     */
};

/*--------------------------------------------------------\
| Name: alnSeqModule [Struct-02:]
| Use:
|  - Structer to hold all function information for alnSeq
\--------------------------------------------------------*/
static struct PyModuleDef alnSeqModule = {
   PyModuleDef_HEAD_INIT,
   "alnSeq", /*Library Name*/
   "alnSeq has a Hirschberg (alnSeqHirsch),\n\
    Needleman Wunsch (alnSeqNeedle), Waterman Smith \n\
    (alnSeqWater), and memory efficent Smith Waterman\n\
    (alnSeqMemWater) pairwise aligners.\n\
      The required input is a reference and query\n\
        sequence. Optional arguments include a\n\
        gap opening score, a gap extension score, and\n\
        a scoring matrix (as file).\n\
      The output is a pair of aligned sequences for\n\
        all algorithims, except the memory efficent\n\
        Smith Waterman, which returns the index 0\n\
        coordinates\n\
   ",
   -1, /*Sets memory to global state. I have no idea here
       ` This is part of subprocess (threading). Right now
       ` -1 sets to not allow. I need to find a good way
       ` to set memory usage.
       */
   alnSeqFunST /*Class with functions calls*/
};

/*--------------------------------------------------------\
| Name: alnSeqModule [Init-01:]
| Use:
|  - Initializes the alnSeq extension for python
\--------------------------------------------------------*/
PyMODINIT_FUNC PyInit_alnSeq(void)
   {return PyModule_Create(&alnSeqModule);}

 /*
   | = arguments after are optional
   $ = comes after |, says that arguments only by keyword
   s = string
   b = uchar (there is no signed character)
   B = uchar
   h = short (Closes you can get to a signed char)
   H = ushort
   i = int
   I = uint
   l = long
   k = ulong (?)
   f = float
   d = double
   p = boolean
*/

