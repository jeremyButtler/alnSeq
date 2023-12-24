/*########################################################
# Name alnSetStruct
# Use:
#  o Holds the stettings structure and functions that
#    support the settings structure
# Includes:
#   - "alnSeqDefaults.h"    (No .c file)
#   - "base10StrToNum.h"    (No .c file)
#   - "dataTypeShortHand.h" (No .c file)
# C Standard libraries:
#   - <stdlib.h>
#   - <stdio.h>
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - Has includes and definede variables
'  o st-01 alnSet:
'     o Holds settings for my alignment program
'  o fun-01 setBpScore:
'    - Sets the score for a base pair (reference/query)
'  o fun-02 setIfBpMatch:
'    - Sets if two bases are a match or mismtach
'  o fun-03 freeAlnSetStack:
'    o Frees variables inside alnSet
'  o fun-04 freeAlnSet:
'    o Frees an alnSet structure (and sets to 0)
'  o fun-05 getBaseScore:
'    - Get the score for a pair of bases from an alignment
'  o fun-06 matchOrSnp:
'    - Check if two bases were a match or mismatch
'  o fun-07 readInScoreFile
'     - Reads in a file of scores for a scoring matrix
'  o fun-08 readInMatchFile:
'    - Reads in a file of matches
'  o fun-09 seqToLookupIndex:
'    - Converts a sequence to a look up table index
'      (table is in alnSetStruct.h)
'  o fun-10 lookupIndexToSeq:
'    - Converts a sequence of lookup indexs back into
'      uppercase characters (a-z)
'  o fun-11 initAlnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Has includes and definede variables
\-------------------------------------------------------*/

#ifndef ALIGNMENT_SETTINGS_STRUCT_H
#define ALIGNMENT_SETTINGS_STRUCT_H

#include <stdlib.h>
#include <stdio.h>

#include "alnSeqDefaults.h"
#include "base10StrToNum.h"
#include "dataTypeShortHand.h"

/*#include <stdio.h>*/

#define defClearNonAlph (1 | 2 | 4 | 8 | 16)
  /*clear 64 bit and case*/
#define defToUper (1 | 2 | 4 | 8 | 16 | 64)
    /*Clear 32nd bit (marks lower case)*/

/*Aligment matrix variables*/
#define defMoveStop 0    /*stop (end of alignment)*/
#define defMoveLeft 1    /*Move left (deletion)*/
#define defMoveUp 2      /*Move up (insertion)*/
#define defMoveDiagnol 3 /*Move on a diagnol (snp/match)*/

/*Flags for match/mismatch*/
#define defBaseMatch 1
#define defBaseSnp 0

#if defined WORDS
   #define defMatrixCol 128
   #undef NOSEQCNVT /*I use a full ascci matrix for 
                    ` sentence/files. There is no
                    ` conversion to be done here
                    */
#elif defined NOSEQCNVT
   #define defMatrixCol 26
#else
   #define defMatrixCol 27
#endif

/*-------------------------------------------------------\
| ST-01: alnSet
| Use: Holds settings for my alignment program
\-------------------------------------------------------*/
typedef struct alnSet
{ /*alnSet*/
   /*Line wrap for printing out an alignment*/
   unsigned int lineWrapUS;
   unsigned short lenFileNameUS;
   char pBasePosBl; /*1 Print out base numbers*/
   char pFullAlnBl;
     /*1: Print out the full alignmnet
     ` 0: Print out the aligned region
     */
   char formatFlag;
     /*defExpandCig: is default format (S D I = X)
     ` defEMBOSS: is EMBOSS format (| space)
     ` defClustal: is clustal format (* space)
     */
   char justScoresBl;
     /*1: print coordiantes & scores only (query/ref scan)
     ` 0: Print the alignment
     */

   /*Preference for alignment algorithim used*/
   char useNeedleBl;
   char useWaterBl;
   char useHirschBl;
   char memWaterBl;
   char twoBitBl;     /*1: use two bit arrays; 0 do not*/

   /*Directional priorities (see alnSeqDefualts.h for
   ` options)
   */
   char bestDirC;

   /*General alignment variables*/
   char noGapBl;     /*1: use gapExtendC; 0 do not*/
   char gapOpenC;    /*Penalty for starting an indel*/
   char gapExtendC;  /* Penalty for extending an indel*/

   /*Waterman smith specific variables*/
   char refQueryScanBl;  /*Keep best score for ref/query*/
     /* If set to 1: Recored a best score for each base
     `  in the reference and query in a Waterman alignment
     */
   char filtQryBl; /*1: filter by query position*/
      /*1: Remove alignments with overlaps in the query
      `    positions
      ` 0: Do not apply this filter
      */
   char filtRefBl; /*1: filter by reference position*/
      /*1: Remove alignments with overlaps in reference
      `    positions
      ` 0: Do not apply this filter
      */
   char filtQryRefBl;
      /*1: Remove alignments with overlaps in both the
      `    query and reference positions
      ` 0: Do not apply this filter
      */
   char pAltAlns;  /*1: Print out alternative alignments*/
   long minScoreL; /*Min score to keep alignment*/

   char scoreMatrixC[defMatrixCol][defMatrixCol];
   char matchMatrixBl[defMatrixCol][defMatrixCol];
     /* Size of both non-WORDS matrixies (27 or 26) is due
     ` to wanting a look up table that can handle
     ` anonymous bases.  Most cells will be set to 0. 
     `  How to get score when not using -DWORDS
     `  score =
     `   scoreMatrixC[(uchar) (base1 & defClearNonAlph)-1]
     `               [(uchar) (base2 & defClearNonAlph)-1]
     */
}alnSet;

/*-------------------------------------------------------\
| Fun-01: setBpScore
|  - Sets the score for a base pair (reference/query)
| Input:
|  - qryBase:
|    o query base to set score fore
|  - refBase:
|    o reference base to set score fore
|  - score:
|    o Score to have for this query/reference combination
|      (max is a short)
|  - alnSetSTPtr:
|    o pointer to alnSet (settings) structer with matrix
|      to add the new query reference comparison score to.
| Output:
|  o Modifies:
|    - one score in an snp/match scoring matrix
|  - WORDS (-DWORDS)
|    o This is the full ascii matrix (128x128)
|  - NOSEQCNVT (-DNOSEQCNVT)
|    o This assumes that the sequences will not be
|      converted to indexes
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
#ifdef WORDS
#define setBpScore(qryBase,refBase,score,alnSetSTPtr){\
    (alnSetSTPtr)->scoreMatrixC\
         [(uchar) (qryBase)]\
         [(uchar) (refBase)]\
       = (score);\
}/*setBpScore with a full ascii matrix*/

#elif defined NOSEQCNVT
#define setBpScore(qryBase,refBase,score,alnSetSTPtr){\
    (alnSetSTPtr)->scoreMatrixC\
         [(uchar) ((qryBase) & defClearNonAlph) - 1]\
         [(uchar) ((refBase) & defClearNonAlph) - 1]\
       = (score);\
}/*setBpScore with no sequence conversion*/

#else

/*The user will not be providing converted bases*/
#define setBpScore(qryBase,refBase,score,alnSetSTPtr){\
   (alnSetSTPtr)->scoreMatrixC\
        [(uchar) ((qryBase) & defClearNonAlph)]\
        [(uchar) ((refBase) & defClearNonAlph)]\
      = (score);\
} /*setBpScore, sequences will be lookup indexes*/
#endif

/*-------------------------------------------------------\
| Fun-02: setIfBpMatch
|  - Sets if two bases are a match or mismtach
| Input:
|  - qryBase:
|    o query base to set score fore
|  - refBase:
|    o reference base to set score fore
|  - match:
|    o 1: the query and reference bases are matches
|    o 0: the query and reference bases are mismatches
|  - alnSetSTPtr:
|    o pointer to alnSet (settings) structer with matrix
|      to add the new query reference comparison score to.
| Output:
|  - Modifies:
|    o one match in the matchMatrix
| Variations:
|  - WORDS (-DWORDS)
|    o This is the full ascii matrix (128x128)
|  - NOSEQCNVT (-DNOSEQCNVT)
|    o This assumes that the sequences will not be
|      converted to indexes
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
#ifdef WORDS
#define setIfBpMatch(qryBase,refBase,match,alnSetSTPtr){\
    (alnSetSTPtr)->matchMatrixBl\
         [(uchar) (qryBase)]\
         [(uchar) (refBase)]\
       = (match);\
}/*setBasePairMatch with a full ascii matrix*/

#elif defined NOSEQCNVT
#define setIfBpMatch(qryBase,refBase,match,alnSetSTPtr){\
    (alnSetSTPtr)->matchMatrixBl\
         [(uchar) ((qryBase) & defClearNonAlph) - 1]\
         [(uchar) ((refBase) & defClearNonAlph) - 1]\
       = (match);\
}/*setBasePairMatch with no sequence conversion*/

#else

/*The user will not be providing converted bases*/
#define setIfBpMatch(qryBase,refBase,match,alnSetSTPtr){\
   (alnSetSTPtr)->matchMatrixBl\
        [(uchar) ((qryBase) & defClearNonAlph)]\
        [(uchar) ((refBase) & defClearNonAlph)]\
      = (match);\
} /*setBasePairMatch, sequences will be lookup indexes*/
#endif

/*-------------------------------------------------------\
| Fun-03: freeAlnSetStack
|  - Does a stack free of an alnSet structer
| Input:
|  - alnSetST:
|    o alnSetST to free internal variables in
|    o alnSetST is not a pointer
| Output:
|  - Free:
|    o Nothing, there are no heap allocated variables. This
|      is here in case I decide to have heap allocated
|      variables on day.
\-------------------------------------------------------*/
#define freeAlnSetStack(alnSetST) ()

/*-------------------------------------------------------\
| Fun-04: freeAlnSet
|  - Frees and alnSet (alignment settings) structure
| Input:
|  - alnSetSTPtr:
|    o Pionter to alnSetST to free
| Output:
|  - Free:
|    o alnSetSTPtr
|  - Set:
|    o alnSetSTPtr to 0 (NULL)
\-------------------------------------------------------*/
#define freeAlnSet(alnSetSTPtr){\
   if(alnSetSTPtr != 0)\
   { /*If: this is a structure to free*/\
      freeAlnSetStack(alnSetSTPtr);\
      free(alnSetSTPtr);\
      alnSetSTPtr = 0;\
   } /*If: this is a structure to free*/\
} /*freeAlnSet*/

/*-------------------------------------------------------\
| Fun-05: getBaseScore:
|  - Get the score for a pair of bases from an alignment
| Input:
|  - qryBase:
|    o character with query base to get score for
|  - refBase:
|    o character with reference base to get score for
|  - alnSetPtr:
|    o alnSet (settings) structer pionter with scoring
|      matrix
| Output:
|  - Returns:
|    o score of the input base pair
| Variations:
|  - WORDS (-DWORDS)
|    o This is the full ascii matrix (128x128)
|  - NOSEQCNVT (-DNOSEQCNVT)
|    o This assumes that the sequences will not be
|      converted to indexes
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
#if defined WORDS
   #define getBaseScore(qryBase, refBase, alnSetSTPtr)(\
      (alnSetSTPtr)->scoreMatrixC\
         [(uchar) ((qryBase)]\
         [(uchar) ((refBase)]\
   ) /*getBasePairScore with no sequence conversion*/

#elif defined NOSEQCNVT
   #define getBaseScore(qryBase, refBase, alnSetSTPtr)(\
      (alnSetSTPtr)->scoreMatrixC\
         [(uchar) ((qryBase) & defClearNonAlph) - 1]\
         [(uchar) ((refBase) & defClearNonAlph) - 1]\
   ) /*getBasePairScore with no sequence conversion*/

#else /*Else the base is pre-converted*/
   #define getBaseScore(qryBase, refBase, alnSetSTPtr)(\
     (alnSetSTPtr)->scoreMatrixC\
        [(uchar) (qryBase)]\
        [(uchar) (refBase)]\
   )/*getBasePairScore, with sequence converted to index*/
#endif

/*-------------------------------------------------------\
| Fun-06: matchOrSnp
|  - Check if two bases were a match or mismatch
| Input:
|  - qryBase:
|    o character with query base to compare
|  - refBase:
|    o character with reference base to compare to
|  - alnSetPtr:
|    o alnSet (settings) structer pionter with scoring
|      matrix
| Output:
|  - Returns:
|    o defBaseMatch if the bases matched (same)
|    o defBaseSnp if bases were not a match (different)
| Variations:
|  - WORDS (-DWORDS)
|    o This is the full ascii matrix (128x128)
|  - NOSEQCNVT (-DNOSEQCNVT)
|    o This assumes that the sequences will not be
|      converted to indexes
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
#if defined WORDS
   #define matchOrSnp(qryBase, refBase, alnSetSTPtr)(\
      (alnSetSTPtr)->matchMatrixBl\
         [(uchar) ((qryBase)]\
         [(uchar) ((refBase)]\
   ) /*getBasePairMatch with no sequence conversion*/

#elif defined NOSEQCNVT
   #define matchOrSnp(qryBase, refBase, alnSetSTPtr)(\
      (alnSetSTPtr)->matchMatrixBl\
         [(uchar) ((qryBase) & defClearNonAlph) - 1]\
         [(uchar) ((refBase) & defClearNonAlph) - 1]\
   ) /*getBasePairMatch with no sequence conversion*/

#else /*Else the base is pre-converted*/
   #define matchOrSnp(qryBase, refBase, alnSetSTPtr)(\
     (alnSetSTPtr)->matchMatrixBl\
        [(uchar) (qryBase)]\
        [(uchar) (refBase)]\
   )/*getBasePairMatch, with sequence converted to index*/
#endif

/*-------------------------------------------------------\
| Fun-07: readInScoreFile
|  - Reads in a file of scores for a scoring matrix
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the score matrix to modify
|  - scoreFILE:
|    o File to get scores from
| Output:
|  - Modifies:
|    o Score matrix in alngSetST to hold the scores from
|      the file (scoreFILE)
|    o scoreFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
static unsigned long readInScoreFile(
    struct alnSet *alnSetSTPtr, /*score matrix to change*/
    FILE *scoreFILE  /*File of scores for scoring matrix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: readInScoreFile
   '  o fun-07 sec-01:
   '    - Variable declerations & set up
   '  o fun-07 sec-02:
   '    - Blank the scoring matrix
   '  o fun-07 sec-03:
   '    - Read in line and check if comment
   '  o fun-07 sec-04:
   '    - Convert score & add to matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ushort lenBuffUS = 1024;
   char buffCStr[lenBuffUS];
   char *tmpCStr = 0;
   char scoreC = 0;

   uchar colUC = 0;
   uchar rowUC = 0;

   buffCStr[lenBuffUS - 1] = '\0';
   buffCStr[lenBuffUS - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < defMatrixCol; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
           alnSetSTPtr->scoreMatrixC[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-07 Sec-03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(fgets(buffCStr, 1024, scoreFILE))
   { /*While I have scores to read in*/
       
       if(buffCStr[0] == '/' && buffCStr[1] == '/')
       { /*On a comment, move onto the next line*/
           while(
               buffCStr[lenBuffUS - 2] != '\0' &&
               buffCStr[lenBuffUS - 2] != '\n'
           ) { /*While have more buffer to read in*/
               buffCStr[lenBuffUS - 2] = '\0';
               fgets(buffCStr, 1024, scoreFILE);
           } /*While have more buffer to read in*/

           /*Reset the buffer*/
           buffCStr[lenBuffUS - 2] = '\0';

           continue;
       } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-07 Sec-04:
       ^  - Convert score & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffCStr[0] == '\n')
           continue;                        /*Blank line*/

       if(buffCStr[0] < 64 && buffCStr[2] < 64)
           return ftell(scoreFILE);  /*Invalid character*/
       
       tmpCStr = base10StrToSC(&buffCStr[4], scoreC);

       setBpScore(
         buffCStr[0],
         buffCStr[2],
         scoreC,
         alnSetSTPtr
       ); /*Add the score to the matrix*/

       if(tmpCStr == &buffCStr[3])
           return ftell(scoreFILE);         /*No score*/

       while(
           buffCStr[lenBuffUS - 2] != '\0' &&
           buffCStr[lenBuffUS - 2] != '\n'
       ){ /*While have more buffer to read in*/
           buffCStr[lenBuffUS - 2] = '\0';
           fgets(buffCStr, 1024, scoreFILE);
       } /*While have more buffer to read in*/

       /*Reset the buffer*/
       buffCStr[lenBuffUS - 2] = '\0';
   } /*While I have scores to read in*/

   return 0;
} /*readInScoreFile*/

/*-------------------------------------------------------\
| Fun-08: readInMatchFile
|  - Reads in a file of matches
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the match matrix to modify
|  - matchFILE:
|    o File to get matchs from
| Output:
|  - Modifies:
|    o Match matrix in alngSetST to hold the matchs from
|      the file (matchFILE)
|    o matchFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
static unsigned long readInMatchFile(
    struct alnSet *alnSetSTPtr,
    FILE *matchFILE  /*File of matchs for scoring matrix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: readInMatchFile
   '  o fun-08 sec-01:
   '    - Variable declerations & set up
   '  o fun-08 sec-02:
   '    - Blank the match matrix
   '  o fun-08 sec-03:
   '    - Read in line and check if comment
   '  o fun-08 sec-04:
   '    - Add match/snp (mismatch) to match matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-08 Sec-01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ushort lenBuffUS = 1024;
   char buffCStr[lenBuffUS];

   uchar colUC = 0;
   uchar rowUC = 0;

   buffCStr[lenBuffUS - 1] = '\0';
   buffCStr[lenBuffUS - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-08 Sec-02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < defMatrixCol; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
           alnSetSTPtr->matchMatrixBl[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-08 Sec-03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(fgets(buffCStr, 1024, matchFILE))
   { /*While I have matchs to read in*/
       
       if(buffCStr[0] == '/' && buffCStr[1] == '/')
       { /*On a comment, move onto the next line*/
           while(
               buffCStr[lenBuffUS - 2] != '\0' &&
               buffCStr[lenBuffUS - 2] != '\n'
           ) { /*While have more buffer to read in*/
               buffCStr[lenBuffUS - 2] = '\0';
               fgets(buffCStr, 1024, matchFILE);
           } /*While have more buffer to read in*/

           /*Reset the buffer*/
           buffCStr[lenBuffUS - 2] = '\0';

           continue;
       } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-08 Sec-04:
       ^  - Convert match & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffCStr[0] == '\n')
           continue;                        /*Blank line*/

       if(buffCStr[4] != '1' && buffCStr[4] != '0')
           return ftell(matchFILE);
           /*This error means I do not know if match/snp*/

       setIfBpMatch(
         buffCStr[0],
         buffCStr[2],
         buffCStr[4],
         alnSetSTPtr
       ); /*Add the match to the matrix*/

       while(
           buffCStr[lenBuffUS - 2] != '\0' &&
           buffCStr[lenBuffUS - 2] != '\n'
       ){ /*While have more buffer to read in*/
           buffCStr[lenBuffUS - 2] = '\0';
           fgets(buffCStr, 1024, matchFILE);
       } /*While have more buffer to read in*/

       /*Reset the buffer*/
       buffCStr[lenBuffUS - 2] = '\0';
   } /*While I have matchs to read in*/

   return 0;
} /*readInMatchFile*/

/*-------------------------------------------------------\
| Fun-09: seqToLookupIndex
|  - Converts a sequence to a look up table index
|    (table is in alnSetStruct.c/h)
| Input:
|  - seqStr:
|    o pointer c-string with the sequence to convert
| Output:
|  - Modifies:
|    o seqST->seqCStr to have look up table indexs (1-27, 
|      with null as 0) instead of bases
| Variations:
|  - WORDS (-DWORDS)
|    - No conversion (blank function)
|  - NOSEQCNVT (-DNOSEQCNVT)
|    - No conversion (blank function)
|  - default:
|    - converts sequence in seqStructPtr to indexes
\-------------------------------------------------------*/
#if defined WORDS || defined NOSEQCNVT
   #define seqToLookupIndex(seqStr) {}
   /*Disable this functino*/
#else
   #define seqToLookupIndex(seqStr){\
      char *tmpStr = seqStr;\
      \
      while(*tmpStr != '\0')\
      { /*Loop: convert bases to lookup table values*/\
         *tmpStr &= defClearNonAlph;\
         ++tmpStr;\
      } /*Loop: convert bases to lookup table values*/\
      \
   } /*seqToLookupIndex*/
#endif

/*-------------------------------------------------------\
| Fun-10: lookupIndexToSeq
|  - Converts a sequence of lookup indexs back into
|    uppercase characters (a-z)
| Input:
|  - seqStructPtr:
|    o Pointer to sequence structer with converte sequence
|      to deconvert (lookup index to base)
| Output:
|  - Modifies:
|    o seqST->seqCStr to have bases instead of look up
|      table indexs
| Variations:
|  - WORDS (-DWORDS)
|    - No conversion (blank function)
|  - NOSEQCNVT (-DNOSEQCNVT)
|    - No conversion (blank function)
|  - default:
|    - converts sequence in seqStructPtr from index to a
|      sequence (human readable)
\-------------------------------------------------------*/
#if defined WORDS || defined NOSEQCNVT
   #define lookupIndexToSeq(seqStr) {}
#else
   #define lookupIndexToSeq(seqStr){\
      char *tmpStr = seqStr;\
      \
      while(*tmpStr != '\0')\
      { /*Loop: convert bases to lookup table values*/\
         *tmpStr |= 64;\
         ++tmpStr;\
      } /*Loop: convert bases to lookup table values*/\
   } /*lookupIndexToSeq*/
#endif /*end of file*/

/*-------------------------------------------------------\
| Fun-11: initAlnSet
|  - Set values in altSet (alingment settings) structure
|    to default values
| Input:
|  - alnSetSTPtr:
|    o poineter to an alnSet (settings) structure to
|      initialize
| Output:
|  o Modifies:
|    - alnSetST to have default alignment settings values
\-------------------------------------------------------*/
static void initAlnSet(
    struct alnSet *alnSetST /*Has settings to initialize*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: initAlnSet
   '  - Set values in altSet (alingment settings)
   '    structure to defaults
   '  o fun-11 sec-01:
   '    - Set non-matrix variables
   '  o fun-11 sec-02:
   '    - Initialize scoring matrix
   '  o fun-11 sec-03:
   '    - Initialize match matrix
   '  o fun-11 sec-04:
   '    - Handle special DNA scoring cases
   '  o fun-11 sec-05:
   '    - Handle special DNA match cases
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-11 Sec-01:
   ^  - Set non-matrix variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Variables for my for loop*/
   uchar colUC = 0;
   uchar rowUC = 0;

   /*Set the default output settings*/
   alnSetST->lineWrapUS = defLineWrap;
   alnSetST->lenFileNameUS = 1024;

   alnSetST->pBasePosBl = defPPos;
   alnSetST->pFullAlnBl = defPAln;
   alnSetST->formatFlag = defFormat;
   alnSetST->justScoresBl = defJustScoresBl;

   /*Choose the aligment style*/
   alnSetST->useNeedleBl = defUseNeedle;
   alnSetST->useWaterBl = defUseWater;
   alnSetST->useHirschBl = defUseHirsch;
   alnSetST->memWaterBl = defUseMemWater;
   alnSetST->twoBitBl = defUseTwoBit;

   /*Select direction to keep if everything is equal*/
   alnSetST->bestDirC = defBestDir;

   /*General alignment variables*/
   alnSetST->noGapBl = defNoGapExtend;
   alnSetST->gapOpenC = defGapOpen;
   alnSetST->gapExtendC = defGapExtend;

   /*Query reference scan variables*/
   alnSetST->refQueryScanBl = defQueryRefScan;
   alnSetST->minScoreL = defMinScore;
   alnSetST->filtQryBl = defFilterByQuery;
   alnSetST->filtRefBl = defFilterByRef;
   alnSetST->filtQryRefBl = defFilterByQueryRef;
   alnSetST->pAltAlns = defPAltAln;
  
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-11 Sec-02:
   ^  - Initialize scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #if defined WORDS
      for(colUC = 0; colUC < 26; ++colUC)
      { /*loop for all columns in the comparison matrix*/
          for(rowUC = 0; rowUC < 26; ++rowUC)
          { /*Loop: Fill in the entire matrix*/ 
             if(colUC != rowUC)
               alnSetST->scoreMatrixC[colUC][rowUC] =
                  defWordSNPScore;
             else
               alnSetST->scoreMatrixC[colUC][rowUC] =
                  defWordMatchScore;
          } /*Loop: Fill in the entire matrix*/ 
      } /*loop for all columns in the comparison matrix*/

   #else  /*Doing the DNA/AA sequence*/
      for(colUC = 0; colUC < defMatrixCol; ++colUC)
      { /*loop for all columns in the comparison matrix*/
          for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
              alnSetST->scoreMatrixC[colUC][rowUC] = 0;
              /*Most of these cells will never be used*/
              /*But are needed to build the table*/
      } /*loop for all columns in the comparison matrix*/
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-11 Sec-03:
   ^  - Initialize match matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Both words, DNA, and AA all are the same when both
    ` characters/bases/amino acids are the same.
    */
    for(colUC = 0; colUC < defMatrixCol; ++colUC)
    { /*loop for all columns in the comparison matrix*/
        for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
        { /*Loop: Fill in the entire matrix*/ 
           if(colUC == rowUC)
             alnSetST->matchMatrixBl[colUC][rowUC] =
                defBaseMatch;
           else
             alnSetST->matchMatrixBl[colUC][rowUC] =
                defBaseSnp;
        } /*Loop: Fill in the entire matrix*/ 
    } /*loop for all columns in the comparison matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-11 Sec-04:
   ^  - Handle special DNA scoring cases
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #if !defined WORDS
      /*Set up scores for non-anonmyous base pairs*/
      setBpScore('a', 'a', defAToA, alnSetST);
      setBpScore('a', 't', defAToT, alnSetST);
      setBpScore('a', 'u', defAToT, alnSetST);
      setBpScore('a', 'g', defAToG, alnSetST);
      setBpScore('a', 'c', defAToC, alnSetST);

      setBpScore('t', 'a', defTToA, alnSetST);
      setBpScore('t', 't', defTToT, alnSetST);
      setBpScore('t', 'g', defTToG, alnSetST);
      setBpScore('t', 'c', defTToC, alnSetST);

      setBpScore('u', 'a', defTToA, alnSetST);
      setBpScore('u', 'u', defTToT, alnSetST);
      setBpScore('u', 'g', defTToG, alnSetST);
      setBpScore('u', 'c', defTToC, alnSetST);

      setBpScore('g', 'a', defGToA, alnSetST);
      setBpScore('g', 't', defGToT, alnSetST);
      setBpScore('g', 'u', defGToT, alnSetST);
      setBpScore('g', 'g', defGToG, alnSetST);
      setBpScore('g', 'c', defGToC, alnSetST);

      setBpScore('c', 'a', defCToA, alnSetST);
      setBpScore('c', 't', defCToT, alnSetST);
      setBpScore('c', 'u', defCToT, alnSetST);
      setBpScore('c', 'g', defCToG, alnSetST);
      setBpScore('c', 'c', defCToC, alnSetST);

      /*non-anonymous base and anonymous base pairs*/
      setBpScore('a', 'w', defAToW, alnSetST);
      setBpScore('a', 's', defAToS, alnSetST);
      setBpScore('a', 'm', defAToM, alnSetST);
      setBpScore('a', 'k', defAToK, alnSetST);
      setBpScore('a', 'r', defAToR, alnSetST);
      setBpScore('a', 'y', defAToY, alnSetST);
      setBpScore('a', 'b', defAToB, alnSetST);
      setBpScore('a', 'd', defAToD, alnSetST);
      setBpScore('a', 'h', defAToH, alnSetST);
      setBpScore('a', 'v', defAToV, alnSetST);
      setBpScore('a', 'n', defAToN, alnSetST);
      setBpScore('a', 'x', defAToX, alnSetST);

      setBpScore('w', 'a', defWToA, alnSetST);
      setBpScore('s', 'a', defSToA, alnSetST);
      setBpScore('m', 'a', defMToA, alnSetST);
      setBpScore('k', 'a', defKToA, alnSetST);
      setBpScore('r', 'a', defRToA, alnSetST);
      setBpScore('y', 'a', defYToA, alnSetST);
      setBpScore('b', 'a', defBToA, alnSetST);
      setBpScore('d', 'a', defDToA, alnSetST);
      setBpScore('h', 'a', defHToA, alnSetST);
      setBpScore('v', 'a', defVToA, alnSetST);
      setBpScore('n', 'a', defNToA, alnSetST);
      setBpScore('x', 'a', defXToA, alnSetST);

      setBpScore('c', 'w', defCToW, alnSetST);
      setBpScore('c', 's', defCToS, alnSetST);
      setBpScore('c', 'm', defCToM, alnSetST);
      setBpScore('c', 'k', defCToK, alnSetST);
      setBpScore('c', 'r', defCToR, alnSetST);
      setBpScore('c', 'y', defCToY, alnSetST);
      setBpScore('c', 'b', defCToB, alnSetST);
      setBpScore('c', 'd', defCToD, alnSetST);
      setBpScore('c', 'h', defCToH, alnSetST);
      setBpScore('c', 'v', defCToV, alnSetST);
      setBpScore('c', 'n', defCToN, alnSetST);
      setBpScore('c', 'x', defCToX, alnSetST);

      setBpScore('w', 'C', defWToC, alnSetST);
      setBpScore('s', 'C', defSToC, alnSetST);
      setBpScore('m', 'C', defMToC, alnSetST);
      setBpScore('k', 'C', defKToC, alnSetST);
      setBpScore('r', 'C', defRToC, alnSetST);
      setBpScore('y', 'C', defYToC, alnSetST);
      setBpScore('b', 'C', defBToC, alnSetST);
      setBpScore('d', 'C', defDToC, alnSetST);
      setBpScore('h', 'C', defHToC, alnSetST);
      setBpScore('v', 'C', defVToC, alnSetST);
      setBpScore('n', 'C', defNToC, alnSetST);
      setBpScore('x', 'C', defXToC, alnSetST);

      setBpScore('g', 'w', defGToW, alnSetST);
      setBpScore('g', 's', defGToS, alnSetST);
      setBpScore('g', 'm', defGToM, alnSetST);
      setBpScore('g', 'k', defGToK, alnSetST);
      setBpScore('g', 'r', defGToR, alnSetST);
      setBpScore('g', 'y', defGToY, alnSetST);
      setBpScore('g', 'b', defGToB, alnSetST);
      setBpScore('g', 'd', defGToD, alnSetST);
      setBpScore('g', 'h', defGToH, alnSetST);
      setBpScore('g', 'v', defGToV, alnSetST);
      setBpScore('g', 'n', defGToN, alnSetST);
      setBpScore('g', 'x', defGToX, alnSetST);

      setBpScore('w', 'g', defWToG, alnSetST);
      setBpScore('s', 'g', defSToG, alnSetST);
      setBpScore('m', 'g', defMToG, alnSetST);
      setBpScore('k', 'g', defKToG, alnSetST);
      setBpScore('r', 'g', defRToG, alnSetST);
      setBpScore('y', 'g', defYToG, alnSetST);
      setBpScore('b', 'g', defBToG, alnSetST);
      setBpScore('d', 'g', defDToG, alnSetST);
      setBpScore('h', 'g', defHToG, alnSetST);
      setBpScore('v', 'g', defVToG, alnSetST);
      setBpScore('n', 'g', defNToG, alnSetST);
      setBpScore('x', 'g', defXToG, alnSetST);

      setBpScore('t', 'w', defTToW, alnSetST);
      setBpScore('t', 's', defTToS, alnSetST);
      setBpScore('t', 'm', defTToM, alnSetST);
      setBpScore('t', 'k', defTToK, alnSetST);
      setBpScore('t', 'r', defTToR, alnSetST);
      setBpScore('t', 'y', defTToY, alnSetST);
      setBpScore('t', 'b', defTToB, alnSetST);
      setBpScore('t', 'd', defTToD, alnSetST);
      setBpScore('t', 'h', defTToH, alnSetST);
      setBpScore('t', 'v', defTToV, alnSetST);
      setBpScore('t', 'n', defTToN, alnSetST);
      setBpScore('t', 'x', defTToX, alnSetST);

      setBpScore('w', 't', defWToT, alnSetST);
      setBpScore('s', 't', defSToT, alnSetST);
      setBpScore('m', 't', defMToT, alnSetST);
      setBpScore('k', 't', defKToT, alnSetST);
      setBpScore('r', 't', defRToT, alnSetST);
      setBpScore('y', 't', defYToT, alnSetST);
      setBpScore('b', 't', defBToT, alnSetST);
      setBpScore('d', 't', defDToT, alnSetST);
      setBpScore('h', 't', defHToT, alnSetST);
      setBpScore('v', 't', defVToT, alnSetST);
      setBpScore('n', 't', defNToT, alnSetST);
      setBpScore('x', 't', defXToT, alnSetST);

      /*Set u & t to same scores (U is RNA version of T)*/
      setBpScore('u', 'w', defTToW, alnSetST);
      setBpScore('u', 's', defTToS, alnSetST);
      setBpScore('u', 'm', defTToM, alnSetST);
      setBpScore('u', 'k', defTToK, alnSetST);
      setBpScore('u', 'r', defTToR, alnSetST);
      setBpScore('u', 'y', defTToY, alnSetST);
      setBpScore('u', 'b', defTToB, alnSetST);
      setBpScore('u', 'd', defTToD, alnSetST);
      setBpScore('u', 'h', defTToH, alnSetST);
      setBpScore('u', 'v', defTToV, alnSetST);
      setBpScore('u', 'n', defTToN, alnSetST);
      setBpScore('u', 'x', defTToX, alnSetST);

      setBpScore('w', 'u', defWToT, alnSetST);
      setBpScore('s', 'u', defSToT, alnSetST);
      setBpScore('m', 'u', defMToT, alnSetST);
      setBpScore('k', 'u', defKToT, alnSetST);
      setBpScore('r', 'u', defRToT, alnSetST);
      setBpScore('y', 'u', defYToT, alnSetST);
      setBpScore('b', 'u', defBToT, alnSetST);
      setBpScore('d', 'u', defDToT, alnSetST);
      setBpScore('h', 'u', defHToT, alnSetST);
      setBpScore('v', 'u', defVToT, alnSetST);
      setBpScore('n', 'u', defNToT, alnSetST);
      setBpScore('x', 'u', defXToT, alnSetST);

      /*anonymous base and anonymous base pairs*/
      setBpScore('w', 'w', defWToW, alnSetST);
      setBpScore('w', 's', defWToS, alnSetST);
      setBpScore('w', 'm', defWToM, alnSetST);
      setBpScore('w', 'k', defWToK, alnSetST);
      setBpScore('w', 'r', defWToR, alnSetST);
      setBpScore('w', 'y', defWToY, alnSetST);
      setBpScore('w', 'b', defWToB, alnSetST);
      setBpScore('w', 'd', defWToD, alnSetST);
      setBpScore('w', 'h', defWToH, alnSetST);
      setBpScore('w', 'v', defWToV, alnSetST);
      setBpScore('w', 'n', defWToN, alnSetST);
      setBpScore('w', 'x', defWToX, alnSetST);

      setBpScore('s', 'w', defSToW, alnSetST);
      setBpScore('s', 's', defSToS, alnSetST);
      setBpScore('s', 'm', defSToM, alnSetST);
      setBpScore('s', 'k', defSToK, alnSetST);
      setBpScore('s', 'r', defSToR, alnSetST);
      setBpScore('s', 'y', defSToY, alnSetST);
      setBpScore('s', 'b', defSToB, alnSetST);
      setBpScore('s', 'd', defSToD, alnSetST);
      setBpScore('s', 'h', defSToH, alnSetST);
      setBpScore('s', 'v', defSToV, alnSetST);
      setBpScore('s', 'n', defSToN, alnSetST);
      setBpScore('s', 'x', defSToX, alnSetST);

      setBpScore('m', 'w', defMToW, alnSetST);
      setBpScore('m', 's', defMToS, alnSetST);
      setBpScore('m', 'm', defMToM, alnSetST);
      setBpScore('m', 'k', defMToK, alnSetST);
      setBpScore('m', 'r', defMToR, alnSetST);
      setBpScore('m', 'y', defMToY, alnSetST);
      setBpScore('m', 'b', defMToB, alnSetST);
      setBpScore('m', 'd', defMToD, alnSetST);
      setBpScore('m', 'h', defMToH, alnSetST);
      setBpScore('m', 'v', defMToV, alnSetST);
      setBpScore('m', 'n', defMToN, alnSetST);
      setBpScore('m', 'x', defMToX, alnSetST);

      setBpScore('k', 'w', defKToW, alnSetST);
      setBpScore('k', 's', defKToS, alnSetST);
      setBpScore('k', 'm', defKToM, alnSetST);
      setBpScore('k', 'k', defKToK, alnSetST);
      setBpScore('k', 'r', defKToR, alnSetST);
      setBpScore('k', 'y', defKToY, alnSetST);
      setBpScore('k', 'b', defKToB, alnSetST);
      setBpScore('k', 'd', defKToD, alnSetST);
      setBpScore('k', 'h', defKToH, alnSetST);
      setBpScore('k', 'v', defKToV, alnSetST);
      setBpScore('k', 'n', defKToN, alnSetST);
      setBpScore('k', 'x', defKToX, alnSetST);

      setBpScore('r', 'w', defRToW, alnSetST);
      setBpScore('r', 's', defRToS, alnSetST);
      setBpScore('r', 'm', defRToM, alnSetST);
      setBpScore('r', 'k', defRToK, alnSetST);
      setBpScore('r', 'r', defRToR, alnSetST);
      setBpScore('r', 'y', defRToY, alnSetST);
      setBpScore('r', 'b', defRToB, alnSetST);
      setBpScore('r', 'd', defRToD, alnSetST);
      setBpScore('r', 'h', defRToH, alnSetST);
      setBpScore('r', 'v', defRToV, alnSetST);
      setBpScore('r', 'n', defRToN, alnSetST);
      setBpScore('r', 'x', defRToX, alnSetST);

      setBpScore('y', 'w', defYToW, alnSetST);
      setBpScore('y', 's', defYToS, alnSetST);
      setBpScore('y', 'm', defYToM, alnSetST);
      setBpScore('y', 'k', defYToK, alnSetST);
      setBpScore('y', 'r', defYToR, alnSetST);
      setBpScore('y', 'y', defYToY, alnSetST);
      setBpScore('y', 'b', defYToB, alnSetST);
      setBpScore('y', 'd', defYToD, alnSetST);
      setBpScore('y', 'h', defYToH, alnSetST);
      setBpScore('y', 'v', defYToV, alnSetST);
      setBpScore('y', 'n', defYToN, alnSetST);
      setBpScore('y', 'x', defYToX, alnSetST);

      setBpScore('b', 'w', defBToW, alnSetST);
      setBpScore('b', 's', defBToS, alnSetST);
      setBpScore('b', 'm', defBToM, alnSetST);
      setBpScore('b', 'k', defBToK, alnSetST);
      setBpScore('b', 'r', defBToR, alnSetST);
      setBpScore('b', 'y', defBToY, alnSetST);
      setBpScore('b', 'b', defBToB, alnSetST);
      setBpScore('b', 'd', defBToD, alnSetST);
      setBpScore('b', 'h', defBToH, alnSetST);
      setBpScore('b', 'v', defBToV, alnSetST);
      setBpScore('b', 'n', defBToN, alnSetST);
      setBpScore('b', 'x', defBToX, alnSetST);

      setBpScore('d', 'w', defDToW, alnSetST);
      setBpScore('d', 's', defDToS, alnSetST);
      setBpScore('d', 'm', defDToM, alnSetST);
      setBpScore('d', 'k', defDToK, alnSetST);
      setBpScore('d', 'r', defDToR, alnSetST);
      setBpScore('d', 'y', defDToY, alnSetST);
      setBpScore('d', 'b', defDToB, alnSetST);
      setBpScore('d', 'd', defDToD, alnSetST);
      setBpScore('d', 'h', defDToH, alnSetST);
      setBpScore('d', 'v', defDToV, alnSetST);
      setBpScore('d', 'n', defDToN, alnSetST);
      setBpScore('d', 'x', defDToX, alnSetST);

      setBpScore('h', 'w', defHToW, alnSetST);
      setBpScore('h', 's', defHToS, alnSetST);
      setBpScore('h', 'm', defHToM, alnSetST);
      setBpScore('h', 'k', defHToK, alnSetST);
      setBpScore('h', 'r', defHToR, alnSetST);
      setBpScore('h', 'y', defHToY, alnSetST);
      setBpScore('h', 'b', defHToB, alnSetST);
      setBpScore('h', 'd', defHToD, alnSetST);
      setBpScore('h', 'h', defHToH, alnSetST);
      setBpScore('h', 'v', defHToV, alnSetST);
      setBpScore('h', 'n', defHToN, alnSetST);
      setBpScore('h', 'x', defHToX, alnSetST);

      setBpScore('v', 'w', defVToW, alnSetST);
      setBpScore('v', 's', defVToS, alnSetST);
      setBpScore('v', 'm', defVToM, alnSetST);
      setBpScore('v', 'k', defVToK, alnSetST);
      setBpScore('v', 'r', defVToR, alnSetST);
      setBpScore('v', 'y', defVToY, alnSetST);
      setBpScore('v', 'b', defVToB, alnSetST);
      setBpScore('v', 'd', defVToD, alnSetST);
      setBpScore('v', 'h', defVToH, alnSetST);
      setBpScore('v', 'v', defVToV, alnSetST);
      setBpScore('v', 'n', defVToN, alnSetST);
      setBpScore('v', 'x', defVToX, alnSetST);

      setBpScore('n', 'w', defNToW, alnSetST);
      setBpScore('n', 's', defNToS, alnSetST);
      setBpScore('n', 'm', defNToM, alnSetST);
      setBpScore('n', 'k', defNToK, alnSetST);
      setBpScore('n', 'r', defNToR, alnSetST);
      setBpScore('n', 'y', defNToY, alnSetST);
      setBpScore('n', 'b', defNToB, alnSetST);
      setBpScore('n', 'd', defNToD, alnSetST);
      setBpScore('n', 'h', defNToH, alnSetST);
      setBpScore('n', 'v', defNToV, alnSetST);
      setBpScore('n', 'n', defNToN, alnSetST);
      setBpScore('n', 'x', defNToX, alnSetST);

      setBpScore('x', 'w', defXToW, alnSetST);
      setBpScore('x', 's', defXToS, alnSetST);
      setBpScore('x', 'm', defXToM, alnSetST);
      setBpScore('x', 'k', defXToK, alnSetST);
      setBpScore('x', 'r', defXToR, alnSetST);
      setBpScore('x', 'y', defXToY, alnSetST);
      setBpScore('x', 'b', defXToB, alnSetST);
      setBpScore('x', 'd', defXToD, alnSetST);
      setBpScore('x', 'h', defXToH, alnSetST);
      setBpScore('x', 'v', defXToV, alnSetST);
      setBpScore('x', 'n', defXToN, alnSetST);
      setBpScore('x', 'x', defXToX, alnSetST);
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-11 Sec-05:
   ^  - Handle special DNA match cases
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #if !defined WORDS
      setIfBpMatch('u', 't', defUMatchT, alnSetST);
      setIfBpMatch('t', 'u', defTMatchU, alnSetST);

      /*non-anonymous base and anonymous base pairs*/
      setIfBpMatch('a', 'w', defAMatchW, alnSetST);
      setIfBpMatch('a', 's', defAMatchS, alnSetST);
      setIfBpMatch('a', 'm', defAMatchM, alnSetST);
      setIfBpMatch('a', 'k', defAMatchK, alnSetST);
      setIfBpMatch('a', 'r', defAMatchR, alnSetST);
      setIfBpMatch('a', 'y', defAMatchY, alnSetST);
      setIfBpMatch('a', 'b', defAMatchB, alnSetST);
      setIfBpMatch('a', 'd', defAMatchD, alnSetST);
      setIfBpMatch('a', 'h', defAMatchH, alnSetST);
      setIfBpMatch('a', 'v', defAMatchV, alnSetST);
      setIfBpMatch('a', 'n', defAMatchN, alnSetST);
      setIfBpMatch('a', 'x', defAMatchX, alnSetST);

      setIfBpMatch('w', 'a', defWMatchA, alnSetST);
      setIfBpMatch('s', 'a', defSMatchA, alnSetST);
      setIfBpMatch('m', 'a', defMMatchA, alnSetST);
      setIfBpMatch('k', 'a', defKMatchA, alnSetST);
      setIfBpMatch('r', 'a', defRMatchA, alnSetST);
      setIfBpMatch('y', 'a', defYMatchA, alnSetST);
      setIfBpMatch('b', 'a', defBMatchA, alnSetST);
      setIfBpMatch('d', 'a', defDMatchA, alnSetST);
      setIfBpMatch('h', 'a', defHMatchA, alnSetST);
      setIfBpMatch('v', 'a', defVMatchA, alnSetST);
      setIfBpMatch('n', 'a', defNMatchA, alnSetST);
      setIfBpMatch('x', 'a', defXMatchA, alnSetST);

      setIfBpMatch('c', 'w', defCMatchW, alnSetST);
      setIfBpMatch('c', 's', defCMatchS, alnSetST);
      setIfBpMatch('c', 'm', defCMatchM, alnSetST);
      setIfBpMatch('c', 'k', defCMatchK, alnSetST);
      setIfBpMatch('c', 'r', defCMatchR, alnSetST);
      setIfBpMatch('c', 'y', defCMatchY, alnSetST);
      setIfBpMatch('c', 'b', defCMatchB, alnSetST);
      setIfBpMatch('c', 'd', defCMatchD, alnSetST);
      setIfBpMatch('c', 'h', defCMatchH, alnSetST);
      setIfBpMatch('c', 'v', defCMatchV, alnSetST);
      setIfBpMatch('c', 'n', defCMatchN, alnSetST);
      setIfBpMatch('c', 'x', defCMatchX, alnSetST);

      setIfBpMatch('w', 'C', defWMatchC, alnSetST);
      setIfBpMatch('s', 'C', defSMatchC, alnSetST);
      setIfBpMatch('m', 'C', defMMatchC, alnSetST);
      setIfBpMatch('k', 'C', defKMatchC, alnSetST);
      setIfBpMatch('r', 'C', defRMatchC, alnSetST);
      setIfBpMatch('y', 'C', defYMatchC, alnSetST);
      setIfBpMatch('b', 'C', defBMatchC, alnSetST);
      setIfBpMatch('d', 'C', defDMatchC, alnSetST);
      setIfBpMatch('h', 'C', defHMatchC, alnSetST);
      setIfBpMatch('v', 'C', defVMatchC, alnSetST);
      setIfBpMatch('n', 'C', defNMatchC, alnSetST);
      setIfBpMatch('x', 'C', defXMatchC, alnSetST);

      setIfBpMatch('g', 'w', defGMatchW, alnSetST);
      setIfBpMatch('g', 's', defGMatchS, alnSetST);
      setIfBpMatch('g', 'm', defGMatchM, alnSetST);
      setIfBpMatch('g', 'k', defGMatchK, alnSetST);
      setIfBpMatch('g', 'r', defGMatchR, alnSetST);
      setIfBpMatch('g', 'y', defGMatchY, alnSetST);
      setIfBpMatch('g', 'b', defGMatchB, alnSetST);
      setIfBpMatch('g', 'd', defGMatchD, alnSetST);
      setIfBpMatch('g', 'h', defGMatchH, alnSetST);
      setIfBpMatch('g', 'v', defGMatchV, alnSetST);
      setIfBpMatch('g', 'n', defGMatchN, alnSetST);
      setIfBpMatch('g', 'x', defGMatchX, alnSetST);

      setIfBpMatch('w', 'g', defWMatchG, alnSetST);
      setIfBpMatch('s', 'g', defSMatchG, alnSetST);
      setIfBpMatch('m', 'g', defMMatchG, alnSetST);
      setIfBpMatch('k', 'g', defKMatchG, alnSetST);
      setIfBpMatch('r', 'g', defRMatchG, alnSetST);
      setIfBpMatch('y', 'g', defYMatchG, alnSetST);
      setIfBpMatch('b', 'g', defBMatchG, alnSetST);
      setIfBpMatch('d', 'g', defDMatchG, alnSetST);
      setIfBpMatch('h', 'g', defHMatchG, alnSetST);
      setIfBpMatch('v', 'g', defVMatchG, alnSetST);
      setIfBpMatch('n', 'g', defNMatchG, alnSetST);
      setIfBpMatch('x', 'g', defXMatchG, alnSetST);

      setIfBpMatch('t', 'w', defTMatchW, alnSetST);
      setIfBpMatch('t', 's', defTMatchS, alnSetST);
      setIfBpMatch('t', 'm', defTMatchM, alnSetST);
      setIfBpMatch('t', 'k', defTMatchK, alnSetST);
      setIfBpMatch('t', 'r', defTMatchR, alnSetST);
      setIfBpMatch('t', 'y', defTMatchY, alnSetST);
      setIfBpMatch('t', 'b', defTMatchB, alnSetST);
      setIfBpMatch('t', 'd', defTMatchD, alnSetST);
      setIfBpMatch('t', 'h', defTMatchH, alnSetST);
      setIfBpMatch('t', 'v', defTMatchV, alnSetST);
      setIfBpMatch('t', 'n', defTMatchN, alnSetST);
      setIfBpMatch('t', 'x', defTMatchX, alnSetST);

      setIfBpMatch('w', 't', defWMatchT, alnSetST);
      setIfBpMatch('s', 't', defSMatchT, alnSetST);
      setIfBpMatch('m', 't', defMMatchT, alnSetST);
      setIfBpMatch('k', 't', defKMatchT, alnSetST);
      setIfBpMatch('r', 't', defRMatchT, alnSetST);
      setIfBpMatch('y', 't', defYMatchT, alnSetST);
      setIfBpMatch('b', 't', defBMatchT, alnSetST);
      setIfBpMatch('d', 't', defDMatchT, alnSetST);
      setIfBpMatch('h', 't', defHMatchT, alnSetST);
      setIfBpMatch('v', 't', defVMatchT, alnSetST);
      setIfBpMatch('n', 't', defNMatchT, alnSetST);
      setIfBpMatch('x', 't', defXMatchT, alnSetST);

      /*Set u & t to same scores (U is RNA version of T)*/
      setIfBpMatch('u', 'w', defTMatchW, alnSetST);
      setIfBpMatch('u', 's', defTMatchS, alnSetST);
      setIfBpMatch('u', 'm', defTMatchM, alnSetST);
      setIfBpMatch('u', 'k', defTMatchK, alnSetST);
      setIfBpMatch('u', 'r', defTMatchR, alnSetST);
      setIfBpMatch('u', 'y', defTMatchY, alnSetST);
      setIfBpMatch('u', 'b', defTMatchB, alnSetST);
      setIfBpMatch('u', 'd', defTMatchD, alnSetST);
      setIfBpMatch('u', 'h', defTMatchH, alnSetST);
      setIfBpMatch('u', 'v', defTMatchV, alnSetST);
      setIfBpMatch('u', 'n', defTMatchN, alnSetST);
      setIfBpMatch('u', 'x', defTMatchX, alnSetST);

      setIfBpMatch('w', 'u', defWMatchT, alnSetST);
      setIfBpMatch('s', 'u', defSMatchT, alnSetST);
      setIfBpMatch('m', 'u', defMMatchT, alnSetST);
      setIfBpMatch('k', 'u', defKMatchT, alnSetST);
      setIfBpMatch('r', 'u', defRMatchT, alnSetST);
      setIfBpMatch('y', 'u', defYMatchT, alnSetST);
      setIfBpMatch('b', 'u', defBMatchT, alnSetST);
      setIfBpMatch('d', 'u', defDMatchT, alnSetST);
      setIfBpMatch('h', 'u', defHMatchT, alnSetST);
      setIfBpMatch('v', 'u', defVMatchT, alnSetST);
      setIfBpMatch('n', 'u', defNMatchT, alnSetST);
      setIfBpMatch('x', 'u', defXMatchT, alnSetST);

      /*anonymous base and anonymous base pairs*/
      setIfBpMatch('w', 's', defWMatchS, alnSetST);
      setIfBpMatch('w', 'm', defWMatchM, alnSetST);
      setIfBpMatch('w', 'k', defWMatchK, alnSetST);
      setIfBpMatch('w', 'r', defWMatchR, alnSetST);
      setIfBpMatch('w', 'y', defWMatchY, alnSetST);
      setIfBpMatch('w', 'b', defWMatchB, alnSetST);
      setIfBpMatch('w', 'd', defWMatchD, alnSetST);
      setIfBpMatch('w', 'h', defWMatchH, alnSetST);
      setIfBpMatch('w', 'v', defWMatchV, alnSetST);
      setIfBpMatch('w', 'n', defWMatchN, alnSetST);
      setIfBpMatch('w', 'x', defWMatchX, alnSetST);

      setIfBpMatch('s', 'w', defSMatchW, alnSetST);
      setIfBpMatch('s', 'm', defSMatchM, alnSetST);
      setIfBpMatch('s', 'k', defSMatchK, alnSetST);
      setIfBpMatch('s', 'r', defSMatchR, alnSetST);
      setIfBpMatch('s', 'y', defSMatchY, alnSetST);
      setIfBpMatch('s', 'b', defSMatchB, alnSetST);
      setIfBpMatch('s', 'd', defSMatchD, alnSetST);
      setIfBpMatch('s', 'h', defSMatchH, alnSetST);
      setIfBpMatch('s', 'v', defSMatchV, alnSetST);
      setIfBpMatch('s', 'n', defSMatchN, alnSetST);
      setIfBpMatch('s', 'x', defSMatchX, alnSetST);

      setIfBpMatch('m', 'w', defMMatchW, alnSetST);
      setIfBpMatch('m', 's', defMMatchS, alnSetST);
      setIfBpMatch('m', 'k', defMMatchK, alnSetST);
      setIfBpMatch('m', 'r', defMMatchR, alnSetST);
      setIfBpMatch('m', 'y', defMMatchY, alnSetST);
      setIfBpMatch('m', 'b', defMMatchB, alnSetST);
      setIfBpMatch('m', 'd', defMMatchD, alnSetST);
      setIfBpMatch('m', 'h', defMMatchH, alnSetST);
      setIfBpMatch('m', 'v', defMMatchV, alnSetST);
      setIfBpMatch('m', 'n', defMMatchN, alnSetST);
      setIfBpMatch('m', 'x', defMMatchX, alnSetST);

      setIfBpMatch('k', 'w', defKMatchW, alnSetST);
      setIfBpMatch('k', 's', defKMatchS, alnSetST);
      setIfBpMatch('k', 'm', defKMatchM, alnSetST);
      setIfBpMatch('k', 'r', defKMatchR, alnSetST);
      setIfBpMatch('k', 'y', defKMatchY, alnSetST);
      setIfBpMatch('k', 'b', defKMatchB, alnSetST);
      setIfBpMatch('k', 'd', defKMatchD, alnSetST);
      setIfBpMatch('k', 'h', defKMatchH, alnSetST);
      setIfBpMatch('k', 'v', defKMatchV, alnSetST);
      setIfBpMatch('k', 'n', defKMatchN, alnSetST);
      setIfBpMatch('k', 'x', defKMatchX, alnSetST);

      setIfBpMatch('r', 'w', defRMatchW, alnSetST);
      setIfBpMatch('r', 's', defRMatchS, alnSetST);
      setIfBpMatch('r', 'm', defRMatchM, alnSetST);
      setIfBpMatch('r', 'k', defRMatchK, alnSetST);
      setIfBpMatch('r', 'y', defRMatchY, alnSetST);
      setIfBpMatch('r', 'b', defRMatchB, alnSetST);
      setIfBpMatch('r', 'd', defRMatchD, alnSetST);
      setIfBpMatch('r', 'h', defRMatchH, alnSetST);
      setIfBpMatch('r', 'v', defRMatchV, alnSetST);
      setIfBpMatch('r', 'n', defRMatchN, alnSetST);
      setIfBpMatch('r', 'x', defRMatchX, alnSetST);

      setIfBpMatch('y', 'w', defYMatchW, alnSetST);
      setIfBpMatch('y', 's', defYMatchS, alnSetST);
      setIfBpMatch('y', 'm', defYMatchM, alnSetST);
      setIfBpMatch('y', 'k', defYMatchK, alnSetST);
      setIfBpMatch('y', 'r', defYMatchR, alnSetST);
      setIfBpMatch('y', 'b', defYMatchB, alnSetST);
      setIfBpMatch('y', 'd', defYMatchD, alnSetST);
      setIfBpMatch('y', 'h', defYMatchH, alnSetST);
      setIfBpMatch('y', 'v', defYMatchV, alnSetST);
      setIfBpMatch('y', 'n', defYMatchN, alnSetST);
      setIfBpMatch('y', 'x', defYMatchX, alnSetST);

      setIfBpMatch('b', 'w', defBMatchW, alnSetST);
      setIfBpMatch('b', 's', defBMatchS, alnSetST);
      setIfBpMatch('b', 'm', defBMatchM, alnSetST);
      setIfBpMatch('b', 'k', defBMatchK, alnSetST);
      setIfBpMatch('b', 'r', defBMatchR, alnSetST);
      setIfBpMatch('b', 'y', defBMatchY, alnSetST);
      setIfBpMatch('b', 'd', defBMatchD, alnSetST);
      setIfBpMatch('b', 'h', defBMatchH, alnSetST);
      setIfBpMatch('b', 'v', defBMatchV, alnSetST);
      setIfBpMatch('b', 'n', defBMatchN, alnSetST);
      setIfBpMatch('b', 'x', defBMatchX, alnSetST);

      setIfBpMatch('d', 'w', defDMatchW, alnSetST);
      setIfBpMatch('d', 's', defDMatchS, alnSetST);
      setIfBpMatch('d', 'm', defDMatchM, alnSetST);
      setIfBpMatch('d', 'k', defDMatchK, alnSetST);
      setIfBpMatch('d', 'r', defDMatchR, alnSetST);
      setIfBpMatch('d', 'y', defDMatchY, alnSetST);
      setIfBpMatch('d', 'b', defDMatchB, alnSetST);
      setIfBpMatch('d', 'h', defDMatchH, alnSetST);
      setIfBpMatch('d', 'v', defDMatchV, alnSetST);
      setIfBpMatch('d', 'n', defDMatchN, alnSetST);
      setIfBpMatch('d', 'x', defDMatchX, alnSetST);

      setIfBpMatch('h', 'w', defHMatchW, alnSetST);
      setIfBpMatch('h', 's', defHMatchS, alnSetST);
      setIfBpMatch('h', 'm', defHMatchM, alnSetST);
      setIfBpMatch('h', 'k', defHMatchK, alnSetST);
      setIfBpMatch('h', 'r', defHMatchR, alnSetST);
      setIfBpMatch('h', 'y', defHMatchY, alnSetST);
      setIfBpMatch('h', 'b', defHMatchB, alnSetST);
      setIfBpMatch('h', 'd', defHMatchD, alnSetST);
      setIfBpMatch('h', 'v', defHMatchV, alnSetST);
      setIfBpMatch('h', 'n', defHMatchN, alnSetST);
      setIfBpMatch('h', 'x', defHMatchX, alnSetST);

      setIfBpMatch('v', 'w', defVMatchW, alnSetST);
      setIfBpMatch('v', 's', defVMatchS, alnSetST);
      setIfBpMatch('v', 'm', defVMatchM, alnSetST);
      setIfBpMatch('v', 'k', defVMatchK, alnSetST);
      setIfBpMatch('v', 'r', defVMatchR, alnSetST);
      setIfBpMatch('v', 'y', defVMatchY, alnSetST);
      setIfBpMatch('v', 'b', defVMatchB, alnSetST);
      setIfBpMatch('v', 'd', defVMatchD, alnSetST);
      setIfBpMatch('v', 'h', defVMatchH, alnSetST);
      setIfBpMatch('v', 'n', defVMatchN, alnSetST);
      setIfBpMatch('v', 'x', defVMatchX, alnSetST);

      setIfBpMatch('n', 'w', defNMatchW, alnSetST);
      setIfBpMatch('n', 's', defNMatchS, alnSetST);
      setIfBpMatch('n', 'm', defNMatchM, alnSetST);
      setIfBpMatch('n', 'k', defNMatchK, alnSetST);
      setIfBpMatch('n', 'r', defNMatchR, alnSetST);
      setIfBpMatch('n', 'y', defNMatchY, alnSetST);
      setIfBpMatch('n', 'b', defNMatchB, alnSetST);
      setIfBpMatch('n', 'd', defNMatchD, alnSetST);
      setIfBpMatch('n', 'h', defNMatchH, alnSetST);
      setIfBpMatch('n', 'v', defNMatchV, alnSetST);
      setIfBpMatch('n', 'x', defNMatchX, alnSetST);

      setIfBpMatch('x', 'w', defXMatchW, alnSetST);
      setIfBpMatch('x', 's', defXMatchS, alnSetST);
      setIfBpMatch('x', 'm', defXMatchM, alnSetST);
      setIfBpMatch('x', 'k', defXMatchK, alnSetST);
      setIfBpMatch('x', 'r', defXMatchR, alnSetST);
      setIfBpMatch('x', 'y', defXMatchY, alnSetST);
      setIfBpMatch('x', 'b', defXMatchB, alnSetST);
      setIfBpMatch('x', 'd', defXMatchD, alnSetST);
      setIfBpMatch('x', 'h', defXMatchH, alnSetST);
      setIfBpMatch('x', 'v', defXMatchV, alnSetST);
      setIfBpMatch('x', 'n', defXMatchN, alnSetST);
   #endif

   return;
} /*initAlnSet*/

#endif
