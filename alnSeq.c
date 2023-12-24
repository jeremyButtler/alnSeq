/*#########################################################
# Name: alignSeq
# Use:
#  - Runs a Needlman Wunsch or Smith Waterman alignment on
#    a pair of fasta files
# Includes:
#  - "hirschberg/hirschberg.h"
#  - "hirschberg/hirschbergNoGap.h"
#  o "hirschberg/genScoreHirsch.h"
#  o "hirschberg/genScoreNoGapHirsc.h"
#  o "hirschberg/genHirsch.h"
#
#  - "needleman/needleman.h"
#  - "needleman/needleNoGap.h"
#  - "needleman/needleTwoBit.h"
#  - "needleman/needleTwoBitNoGap.h"
#  o "needleman/genNeedle.h"
#  o "needleman/genNeedleNoGap.h"
#
#  - "waterman/waterman.h"
#  - "waterman/waterNoGap.h"
#  - "waterman/waterTwoBit.h"
#  - "waterman/waterTwoBitNoGap.h"
#  o "waterman/genWater.h"
#  o "waterman/genWaterNoGap.h"
#
#  - "waterman/waterScan.h"
#  - "waterman/waterScanNoGap.h"
#  - "waterman/waterScanTwoBit.h"
#  - "waterman/waterScanTwoBitNoGap.h"
#  o "waterman/genWaterScan.h"
#  o "waterman/genWaterScanNoGap.h"
#
#  - "memWater/memWater.h"
#  - "memWater/memWaterNoGap.h"
#
#  - "memWater/memWaterScan.h"
#  - "memWater/memWaterScanNoGap.h"
#
#  - "general/sortAndFiltAltAlns.h"
#  o "general/alnMatrixStruct.h"
#  o "general/alnSeqDefaults.h"
#  o "general/alnSetStruct.h"
#  o "general/alnStruct.h"
#  o "general/base10StrToNum.h"
#  o "general/dataTypeShortHand.h"
#  o "general/genAln.h"
#  o "general/genScan.h"
#  o "general/seqStruct.h"
#  o "general/twoBitArrays.h"
#  o "general/genMath.h"
# C standard libraries:
#  o <string.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <stdint.h>
#########################################################*/

#include "hirschberg/hirschberg.h"
#include "hirschberg/hirschbergNoGap.h"

#include "memWater/memWater.h"
#include "memWater/memWaterNoGap.h"

#include "memWater/memWaterScan.h"
#include "memWater/memWaterScanNoGap.h"

#include "waterman/waterman.h"
#include "waterman/watermanNoGap.h"
#include "waterman/waterTwoBit.h"
#include "waterman/waterTwoBitNoGap.h"

#include "waterman/waterScan.h"
#include "waterman/waterScanNoGap.h"
#include "waterman/waterScanTwoBit.h"
#include "waterman/waterScanTwoBitNoGap.h"

#include "needleman/needleman.h"
#include "needleman/needleNoGap.h"
#include "needleman/needleTwoBit.h"
#include "needleman/needleTwoBitNoGap.h"

#include "general/sortAndFiltAltAlns.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOP: Start Of Program
'  - main: Run this program
'  o fun-01 checkInput:
'    o Gets user input
'  o fun-02 printHelpMesg:
'    - Prints the help message for alnSeq
'  o fun-03 printCompilerSettings:
'    - Prints out the compiler flags used
'    - Prints out what each possible compiler flag does
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Fun-01: checkInput
|  - Gets user input and lets program now if something
|    was incorrect.
| Input:
|  - lenArgsInt:
|    o Number of arguments the user input
|  - argsCStr:
|    o Array of c-strings with the users input arguments
|  - refFileCStr:
|    o Pointer to c-string to hold the name of the fasta
|      file with the reference sequence
|  - qryFileCStr:
|    o Pointer to c-string to hold the name of the fasta
|      file with the query sequence
|  - outFileStr:
|    o Pointer to c-string to hold the name of the file to
|      output the main alignment to
|  - altAlnFileStr:
|    o Pointer to c-string to hold the name of the file to
|      output the alternative alignments to
|  - refStartAlnUL
|    o Pointer to unsigned long to hold the first position
|      to align on the reference sequence
|  - refEndAlnUL
|    o Pointer to unsigned long to hold the last position
|      to align on the reference sequence
|  - qryStartAlnUL
|    o Pointer to unsigned long to hold the first position
|      to align on the query sequence
|  - qryEndAlnUL
|    o Pointer to unsigned long to hold the last position
|      to align on the query sequence
|  - scoreMtrxFileStr:
|    o Pointer to hold the file name of the input score
|      matrix (socre matrix is read into settings by
|      checkInput)
|  - matchMtrxFileStr:
|    o Pointer to hold the file name of the input match
|      matrix (the match matrix is read into settings by
|      checkInput)
|  - settings:
|    o Pointer to alnSet structure to hold the users
|      input settings for the alignment
| Output:
|  - Modifies:
|    o Each input variable to hold the specified user
|      input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to "-score-matrix" for invalid scoring
|       matrix input
|  - Prints to stdout when the scoring file is invalid
\-------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,     /*Number arguments user input*/
    char *argsCStr[],    /*Array with user arguments*/
    char **refFileCStr,  /*file name of reference file*/
    char **queryFileCStr,/*File name of the query file*/
    char **outFileStr,     /*Name of the output file*/
    char **altAlnFileStr, /*file: alternative alignments*/
    ulong *refStartAlnUL, /*Start of reference alignment*/
    ulong *refEndAlnUL,   /*End of reference alignment*/
    ulong *qryStartAlnUL, /*Start of query alignment*/
    ulong *qryEndAlnUL,   /*End of query alignment*/
    char **scoreMtrxFileStr,/*Holds scoring matrix file*/
    char **matchMtrxFileStr,/*Holds matching matrix file*/
    struct alnSet *settings /*Aligment settings*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '  o fun-01 sec-01:
   '    -  Variable declerations
   '  o fun-01 sec-02:
   '    - Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Fun-02: printHelpMesg
|  - Prints the help message for alnSeq
| Input:
|  - outFILE:
|    o File to print the help message to
|  - breifBl:
|    o 1: Print the short help message
|    o 0: Print the entire help message
| Output:
|  - Prints:
|    o help message to outFILE.
\-------------------------------------------------------*/
void printHelpMesg(
   FILE *outFILE,
   char breifBl   /*Print a shorter help message*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: printHelpMesg
   '  - Prints out the help message to outFILE
   '  o fun-02 sec-01:
   '    - Usage block
   '  o fun-02 sec-02:
   '    - Input block
   '  o fun-02 sec-03:
   '    - Output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Fun-03: printCompilerSettings
|   - Prints out the compiler flags used
|   - Prints out what each possible compiler flag does
| Input:
|   - outFILE:
|     o File to print flags and flag descriptions to
|   - pDescBl:
|     o 1: print out descriptions for each flag (all)
|     o 0: Do not print out any flag descriptions
| Output:
|   - Prints flags and description to outFILE
\-------------------------------------------------------*/
void printCompileSettings(
   FILE *outFILE, /*Output file*/
   char pDescBl   /*not 0: print out flag descriptions*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC:
   '  - Print out compile flags used and what each
   '    compiler flag does.
   '  o fun-03 sec-01:
   '    - Print out the compiled settings
   '  o fun-03 sec-02:
   '    - Print out what each compiler flag does
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

int main(
    int lenArgsInt,
    char *argsCStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '  - Main function that drives the Smith Waterman or
   '    Needlman Wunsch alignment
   '  o main sec-01:
   '    - Variable declerations
   '  o main sec-02:
   '    - Read in user input and check input
   '  o main sec-03:
   '    - read in the reference sequence
   '  o main sec-04:
   '    - read in the query sequence
   '  o main sec-05:
   '    - Do the alingment
   '  o main sec-06:
   '    - Print out multi-alignments or find the single
   '      alignment array
   '  o main sec-07:
   '    - Print out the alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*User input*/
   char *refFileCStr = 0;
   char *queryFileCStr = 0;
   char *outFileStr = 0;
   char *altAlnFileStr = 0;
   char *inputCStr = 0;
   char *scoreMtrxFileStr = 0;
   char *matchMtrxFileStr = 0;

   /*Holds the reference and query sequence*/
   struct seqStruct refST;
   struct seqStruct queryST;

   /*Caputures error type from functions*/
   uchar errUC = 0;
   long bestScoreL = 0;
   ulong iterUL = 0;

   /*Will hold the users settings*/
   struct alnSet settings;

   /*For holding alignment output*/
   struct alnMatrix *alnMtrxST = 0;
   struct alnMatrixTwoBit *alnMtrxTwoBitST = 0;
   struct alnStruct *alnST = 0;

   FILE *faFILE = 0;
   FILE *outFILE = 0;
   FILE *altAlnFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-02:
   ^  - initalize variables, get user input, & check input
   ^  o main sec-02 sub-01:
   ^    - initialize variables and get user input
   ^  o main sec-02 sub-02:
   ^    - Check user input for errors
   ^  o main sec-02 sub-03:
   ^    - Check if the main output file exists/was input
   ^  o main sec-02 sub-04:
   ^    - Check if file input for alternative alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec-02 Sub-01:
   *  - Initialize variables and Get user input
   \*****************************************************/

   initSeqST(&refST);
   initSeqST(&queryST);
   initAlnSet(&settings);

   inputCStr =
       checkInput(
           &lenArgsInt,
           argsCStr,
           &refFileCStr,
           &queryFileCStr,
           &outFileStr,
           &altAlnFileStr, /*TODO: ADD IN SETTING*/
           &(refST.offsetUL), /*TODO: ADD IN CALL*/
           &(refST.endAlnUL),
           &(queryST.offsetUL),
           &(queryST.endAlnUL),
           &scoreMtrxFileStr,
           &matchMtrxFileStr,
           &settings
   ); /*Get user input*/

   /*****************************************************\
   * Main Sec-02 Sub-02:
   *  - Check user input for errors
   \*****************************************************/

   if(inputCStr != 0)
   { /*If had problematic input*/
        if(strcmp(inputCStr, "-h") == 0 ||
           strcmp(inputCStr, "--h") == 0 ||
           strcmp(inputCStr, "-help") == 0 ||
           strcmp(inputCStr, "--help") == 0 ||
           strcmp(inputCStr, "help") == 0
        ) { /*If user wanted the help message*/
            printHelpMesg(stdout, 1);/*Breif help messag*/
            exit(0);
        } /*If user wanted the help message*/

        if(strcmp(inputCStr, "-h-all") == 0 ||
           strcmp(inputCStr, "--h-all") == 0 ||
           strcmp(inputCStr, "-help-all") == 0 ||
           strcmp(inputCStr, "--help-all") == 0 ||
           strcmp(inputCStr, "help-all") == 0
        ) { /*If user wanted the help message*/
            printHelpMesg(stdout, 0); /*full help messag*/
            exit(0);
        } /*If user wanted the help message*/


        if(strcmp(inputCStr, "-V") == 0 ||
           strcmp(inputCStr, "-v") == 0 ||
           strcmp(inputCStr, "--V") == 0 ||
           strcmp(inputCStr, "--v") == 0 ||
           strcmp(inputCStr, "--version") == 0 ||
           strcmp(inputCStr, "--Version") == 0 ||
           strcmp(inputCStr, "-version") == 0 ||
           strcmp(inputCStr, "-Version") == 0 ||
           strcmp(inputCStr, "Version") == 0 ||
           strcmp(inputCStr, "version") == 0
        ) { /*if the user wanted the version number*/
            fprintf(
                stdout,
                "alnSeq version: %u\n",
                defVersion
            ); /*Print out closest thing to a version*/
            exit(0);
        } /*Else if the user wanted the version number*/

        else if(strcmp(inputCStr, "-flags") == 0)
        { /*If user wanted to print the compiler flags*/
           printCompileSettings(stdout, 1);
           exit(0);
        } /*If user wanted to print the compiler flags*/

        /*Just the flags (no descriptions)*/
        else if(strcmp(inputCStr, "-flags-only") == 0)
        { /*If user wanted to print the compiler flags*/
            fprintf(
                stdout,
                "alnSeq version: %u\n",
                defVersion
            ); /*user will likely want version number*/

           printCompileSettings(stdout, 0);
           exit(0);
        } /*If user wanted to print the compiler flags*/

        /*If file with scoring matrix was invalid*/
        else if(strcmp(inputCStr, "-score-matrix") == 0)
            exit(1);

        else if(strcmp(inputCStr, "-match-matrix") == 0)
            exit(1);

        else if(inputCStr != 0)
        { /*If user had invalid input*/
            printHelpMesg(stderr, 1); /*short help*/
            fprintf(stderr, "%s is invalid\n", inputCStr);
            exit(1); /*Let user know their was an error*/
        } /*If user had invalid input*/
   } /*If had problematic input*/;

   /*****************************************************\
   * Main Sec-02 Sub-03:
   *  - Check if the main output file exists/was input
   \*****************************************************/

   if(outFileStr != 0)
   { /*If printing output to a file*/
        outFILE = fopen(outFileStr, "w");

        if(outFILE == 0)
        { /*If an invalid output file*/
            printf(
              "Output (-out %s) file is invalid.\n",
              outFileStr
            ); /*Let user know about the invalid file*/

            exit(-1);
        } /*If an invalid output file*/

        fclose(outFILE);
        outFILE = 0;
   } /*If printing output to a file*/

   /*****************************************************\
   * Main Sec-02 Sub-04:
   *  - Check if file was input for alternative alignments
   \*****************************************************/

   if(altAlnFileStr != 0)
   { /*If: I have a file for alternative alignments*/
        altAlnFILE = fopen(altAlnFileStr, "w");

        if(altAlnFILE == 0)
        { /*If an invalid output file*/
            printf(
              "Output (-out-alt %s) file is invalid.\n",
              altAlnFileStr
            ); /*Let user know about the invalid file*/

            exit(-1);
        } /*If an invalid output file*/

        fclose(altAlnFILE);
        altAlnFILE = 0;
   } /*If: I have a file for alternative alignments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-3:
   ^  - read in the reference sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   faFILE = fopen(refFileCStr, "r");

   if(faFILE == 0) 
   { /*If: reference file could not be opened*/
       fprintf(
         stderr,
         "Reference (-ref %s) could not be opend\n",
         refFileCStr
       );

       exit(-1);
   } /*If: reference file could not be opened*/

   /*Read in the reference sequence*/
   errUC = readFaSeq(faFILE, &refST);
   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { /*If: I had an Invalid fasta file*/
       freeSeqSTStack(&refST);

       fprintf(
         stderr,
         "Reference (-ref %s) is not valid\n",
         refFileCStr
       );

       exit(-1);
   } /*If: I had an Invalid fasta file*/

   if(errUC & 64)
   { /*If: I had a memory error*/
       freeSeqSTStack(&refST);
       fprintf(stderr, "Memory allocation error\n");
       exit(-1);
   } /*If: I had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-04:
   ^  - read in the query sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   faFILE = fopen(queryFileCStr, "r");

   if(faFILE == 0) 
   { /*If: reference file could not be opened*/
       freeSeqSTStack(&refST);

       fprintf(
         stderr,
         "Query (-query %s) could not be opend\n",
         queryFileCStr
       );

       exit(-1);
   } /*If: reference file could not be opened*/


   /*Read in the query sequence*/
   errUC = readFaSeq(faFILE, &queryST);
   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { /*If: I had an Invalid fasta file*/
       freeSeqSTStack(&refST);
       freeSeqSTStack(&queryST);

       fprintf(
         stderr,
         "Query (-query %s) is not valid\n",
         refFileCStr
       );
       exit(-1);
   } /*If: I had an Invalid fasta file*/

   if(errUC & 64)
   { /*If: I had a memory error*/
      freeSeqSTStack(&refST);
      freeSeqSTStack(&queryST);
      fprintf(stderr, "Memory allocation error\n");
      exit(-1);
   } /*If: I had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-05:
   ^  - Do the alingment
   ^  o main sec-05 sub-01:
   ^    - Set up for the alignment
   ^  o main sec-05 sub-02:
   ^    - Open the output files
   ^  o main sec-05 sub-03:
   ^    - Check if doing a Needleman alignment
   ^  o main sec-05 sub-04:
   ^    - Check if doing an Waterman alignment
   ^  o main sec-05 sub-05:
   ^    - Check if doing an query reference scan Waterman
   ^  o main sec-05 sub-06:
   ^    - Check if doing an memory efficent Waterman
   ^      alignment that prints out alternative alignment
   ^      positions (Gets combined with a Hirschberg)
   ^  o main sec-05 sub-07:
   ^    - Check if doing an memory efficent Waterman
   ^      alignment that prints out a single alignment
   ^      (uses a Hirschberg to get the alignment)
   ^  o main sec-05 sub-08:
   ^    - Check if doing an Hirschberg alignment
   ^  o main sec-05 sub-09:
   ^    - Let user know this was an invalid alignment
   ^  o main sec-05 sub-10:
   ^    - Check if the alignment failed
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Main Sec-05 Sub-01:
   *  - Set up for the alignment
   \******************************************************/

   seqToLookupIndex(refST.seqCStr);
   seqToLookupIndex(queryST.seqCStr);

   if(refST.endAlnUL == 0)
      refST.endAlnUL = refST.lenSeqUL - 1;
   if(queryST.endAlnUL == 0)
      queryST.endAlnUL = queryST.lenSeqUL - 1;

   /*****************************************************\
   *  Main Sec-05 Sub-02:
   *   - Open the output files
   \*****************************************************/

   if(outFileStr != 0)
   { /*If: The user input an alternative alignment file*/
      if(strcmp("-", outFileStr) == 0)
         outFILE = stdout;
      else outFILE = fopen(outFileStr, "w");
   } /*If: The user input an alternative alignment file*/

   else
   { /*Else putting output to stdout*/
      outFILE = stdout;
      outFileStr = "out";
   } /*Else putting output to stdout*/

   if(altAlnFileStr != 0)
   { /*If: The user input an alternative alignment file*/
      if(strcmp("-", altAlnFileStr) == 0)
         altAlnFILE = stdout;
      else altAlnFILE = fopen(altAlnFileStr, "w");
   } /*If: The user input an alternative alignment file*/

   else altAlnFILE = outFILE;

   /******************************************************\
   * Main Sec-05 Sub-03:
   *  - Check if doing an Needleman alignment
   \******************************************************/

   if(settings.useNeedleBl)
   { /*If: I am doing a Needleman alignment*/
     if(settings.noGapBl)
        alnMtrxST =
           NeedleAlnNoGap(&queryST, &refST, &settings);
     else if(settings.twoBitBl && settings.noGapBl)
        alnMtrxTwoBitST = 
           NeedleTwoBitNoGap(&queryST, &refST, &settings);
     else if(settings.twoBitBl)
        alnMtrxTwoBitST = 
           NeedleTwoBit(&queryST, &refST, &settings);
     else
        alnMtrxST =
           NeedlemanAln(&queryST, &refST, &settings);

     if(alnMtrxST == 0 && alnMtrxTwoBitST == 0)
     { /*If: the aligment falied*/
         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);

         fprintf(
            stderr,
           "Ran out of memory for Needleman alignment\n"
         );

         if(outFILE != stdout) fclose(outFILE);
         if(altAlnFILE != stdout) fclose(altAlnFILE);

         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);

         exit(-1);
     } /*If: the aligment falied*/

     if(alnMtrxST != 0)
     { /*If: I am doning a byte array (normal)*/
         alnST =
            dirMatrixToAln(
               &refST,
               &queryST,
               alnMtrxST->bestEndIndexUL,
               &settings,
               alnMtrxST
         );

        bestScoreL = alnMtrxST->bestScoreL;
        freeAlnMatrix(alnMtrxST); /*No longer needed*/
        alnMtrxST = 0;
     } /*If: I am doning a byte array (normal)*/

     else
     { /*Else: I am doing a two bit alignment*/
         alnST =
            twoBitDirMatrixToAln(
               &refST,
               &queryST,
               alnMtrxTwoBitST->bestEndIndexUL,
               &settings,
               alnMtrxTwoBitST
         );

        bestScoreL = alnMtrxTwoBitST->bestScoreL;
        freeAlnMatrixTwoBit(alnMtrxTwoBitST);
        alnMtrxTwoBitST = 0;
     } /*Else: I am doing a two bit alignment*/

     goto printAlignment;
   } /*If: I am doing a Needleman alignment*/

   /******************************************************\
   * Main Sec-05 Sub-03:
   *  - Check if doing an Waterman alignment
   \******************************************************/

   else if(settings.useWaterBl && !settings.refQueryScanBl)
   { /*Else If: doing a waterman alignment*/
      if(settings.noGapBl)
         alnMtrxST =
            WatermanAlnNoGap(&queryST, &refST, &settings);
      else if(settings.twoBitBl && settings.noGapBl)
         alnMtrxTwoBitST =
            WaterTwoBitNoGap(&queryST, &refST, &settings);
      else if(settings.twoBitBl)
         alnMtrxTwoBitST =
            WaterTwoBit(&queryST, &refST, &settings);
      else
         alnMtrxST =
            WatermanAln(&queryST,&refST,&settings);

      if(alnMtrxST == 0 && alnMtrxTwoBitST == 0)
      { /*If: the aligment falied*/
          freeSeqSTStack(&refST);
          freeSeqSTStack(&queryST);
  
         fprintf(
             stderr,
             "Ran out of memory for Waterman alignment\n"
          );
  
         if(outFILE != stdout) fclose(outFILE);
         if(altAlnFILE != stdout) fclose(altAlnFILE);
   
         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);

         exit(-1);
      } /*If: the aligment falied*/

      if(alnMtrxST != 0)
      { /*If: I did a byte Waterman Smith alignment*/
          alnST =
             dirMatrixToAln(
                &refST,
                &queryST,
                alnMtrxST->bestEndIndexUL,
                &settings,
                alnMtrxST
          );

         bestScoreL = alnMtrxST->bestScoreL;
         freeAlnMatrix(alnMtrxST); /*No longer need*/
         alnMtrxST = 0;
      } /*If: I did a byte Waterman Smith alignment*/

      else if(alnMtrxTwoBitST != 0)
      { /*If: I did a byte Waterman Smith alignment*/
          alnST =
             twoBitDirMatrixToAln(
                &refST,
                &queryST,
                alnMtrxTwoBitST->bestEndIndexUL,
                &settings,
                alnMtrxTwoBitST
          );

         bestScoreL = alnMtrxTwoBitST->bestScoreL;
         freeAlnMatrixTwoBit(alnMtrxTwoBitST);
         alnMtrxTwoBitST = 0;
      } /*If: I did a byte Waterman Smith alignment*/

      goto printAlignment;
   } /*Else If: I am doing a single score waterman*/

   /*****************************************************\
   * Main Sec-05 Sub-04:
   *  - Check if doing an query reference scan Waterman
   \*****************************************************/

   else if(settings.useWaterBl && settings.refQueryScanBl)
   { /*Else If: Wateman with query reference scan */
      if(settings.noGapBl)
         alnMtrxST =
            WaterScanNoGap(&queryST, &refST, &settings);
      else if(settings.twoBitBl && settings.noGapBl)
         alnMtrxTwoBitST =
          WaterScanTwoBitNoGap(&queryST,&refST,&settings);
      else if(settings.twoBitBl)
         alnMtrxTwoBitST =
            WaterScanTwoBit(&queryST, &refST, &settings);
      else
         alnMtrxST =
            WaterScan(&queryST, &refST, &settings);

      if(alnMtrxST == 0 && alnMtrxTwoBitST == 0)
      { /*If: the aligment falied*/
         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);
    
         fprintf(
              stderr,
              "Ran out of memory for Waterman Query"
         );
         fprintf(stderr, " reference scan\n");
    
         if(outFILE != stdout) fclose(outFILE);
         if(altAlnFILE != stdout) fclose(altAlnFILE);

         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);

         exit(-1);
      } /*If: the aligment falied*/
  
      /*Check if I am filtering the alternative scores*/
      /*
      if(settings.filtQryRefBl)
         altScoreFiltRefQry(alnMtrxST);
      if(settings.filtRefBl)
         altScoreFiltRef(alnMtrxST);
      if(settings.filtQryBl)
         altScoreFiltQry(alnMtrxST);
      */
      if(!settings.pAltAlns)
      { /*If: I am not printing full alt alignments*/
         if(alnMtrxST != 0)
         { /*If: I used a byte matrix*/
            pAltAlnScores(
               alnMtrxST,
               settings.minScoreL,
               altAlnFILE
            ); /*Print out socres for all alignments*/
         } /*If: I used a byte matrix*/

         else
         { /*Else: I used a two bit matrix*/
            pAltAlnScores(
               alnMtrxTwoBitST,
               settings.minScoreL,
               altAlnFILE
            ); /*Print out socres for all alignments*/
         } /*Else: I used a two bit matrix*/

         if(settings.justScoresBl)
         { /*If: I am printing out scores only*/
            if(outFILE != stdout) fclose(outFILE);
            if(altAlnFILE != stdout) fclose(altAlnFILE);

            freeSeqSTStack(&refST);
            freeSeqSTStack(&queryST);

            if(alnMtrxST != 0)
              {freeAlnMatrix(alnMtrxST);}
            else
              {freeAlnMatrixTwoBit(alnMtrxTwoBitST);}
   
            exit(0);
         } /*If: I am printing out scores only*/
      } /*If: I am not printing full alt alignments*/

      else
      { /*Else: printing out an alignment for each alt*/
         for(
            iterUL = 0;
            iterUL < alnMtrxST->lenArraysUL;
            ++iterUL
         ){ /*Loop: Pint out kept alignments*/

            if(alnMtrxST != 0)
            { /*If: I am doing a byte alignment*/
               alnST =
                  dirMatrixToAln(
                     &refST,
                     &queryST,
                     alnMtrxST->endIndexAryUL[iterUL],
                     &settings,
                     alnMtrxST
               ); /*Get the alternative alignment*/
            } /*If: I am doing a byte alignment*/

            else
            { /*Else: I am doing a two bit alignment*/
               alnST =
                 twoBitDirMatrixToAln(
                   &refST,
                   &queryST,
                   alnMtrxTwoBitST->endIndexAryUL[iterUL],
                   &settings,
                   alnMtrxTwoBitST
               ); /*Get the alternative alignment*/
            } /*Else: I am doing a two bit alignment*/
 
            errUC =             
               printAln(
                  altAlnFILE,
                  altAlnFileStr,
                  &refST,
                  &queryST,
                  alnST,
                  alnMtrxST->scoreAryL[iterUL],
                  &settings,
                  scoreMtrxFileStr/*Name: scoring matrix*/
            ); /*Print out the alternative alignment*/

            freeAlnST(alnST);
            alnST = 0;

            if(errUC)
            { /*If: I falied to print out the alignment*/
               if(altAlnFILE != stdout)
                  fclose(altAlnFILE);

               if(outFILE != stdout) fclose(outFILE);

               freeSeqSTStack(&refST);
               freeSeqSTStack(&queryST);

               if(alnMtrxST != 0)
                 {freeAlnMatrix(alnMtrxST);}
               else
                 {freeAlnMatrixTwoBit(alnMtrxTwoBitST);}
   
               fprintf(
                  stderr,
                  "Error while printing alternative"
               );
               fprintf(stderr,  "alignments\n");

               exit(-1);
            } /*If: I falied to print out the alignment*/
         } /*Loop: Pint out kept alignments*/
      } /*Else: printing out an alignment for each alt*/
  
      if(alnMtrxST != 0)
      { /*If: I did a byte aligment*/
         alnST =
            dirMatrixToAln(
               &refST,
               &queryST,
               alnMtrxST->bestEndIndexUL,
               &settings,
               alnMtrxST
         );
  
         bestScoreL = alnMtrxST->bestScoreL;
         freeAlnMatrix(alnMtrxST); /*No longer need*/
      } /*If: I did a byte aligment*/

      else
      { /*Else: I did a two bit aligment*/
         alnST =
            twoBitDirMatrixToAln(
               &refST,
               &queryST,
               alnMtrxTwoBitST->bestEndIndexUL,
               &settings,
               alnMtrxTwoBitST
         );
  
         bestScoreL = alnMtrxTwoBitST->bestScoreL;
         freeAlnMatrixTwoBit(alnMtrxTwoBitST);
      } /*Else: I did a two bit aligment*/

      goto printAlignment;
   } /*Else If: Wateman with query reference scan*/

   /******************************************************\
   * Main Sec-05 Sub-04:
   *  - Check if doing an memory efficent Waterman
   *    alignment that prints out alternative alignment
   *    positions (Gets combined with a Hirschberg)
   \******************************************************/

   else if(settings.memWaterBl && settings.refQueryScanBl)
   { /*Else if; memory waterman alignment + scan*/
     if(settings.noGapBl)
        alnMtrxST =
           memWaterScanNoGap(&queryST,&refST,&settings);
     else
        alnMtrxST =
           memWaterScan(&queryST,&refST,&settings);

      if(alnMtrxST == 0)
      { /*If: the aligment falied*/
         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);
       
         fprintf(
            stderr,
            "Ran out of memory for mem-water Query"
         );
         fprintf(stderr, " reference scan\n");
       
         if(outFILE != stdout) fclose(outFILE);
         if(altAlnFILE != stdout) fclose(altAlnFILE);
   
         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);
   
         exit(-1);
      } /*If: the aligment falied*/


      /*Check if I am filtering the alternative scores*/
      /*
      if(settings.filtQryRefBl)
         altScoreFiltRefQry(alnMtrxST);
      if(settings.filtRefBl)
         altScoreFiltRef(alnMtrxST);
      if(settings.filtQryBl)
         altScoreFiltQry(alnMtrxST);
      */
       alnMatrixSortScores(alnMtrxST, 0, alnMtrxST->lenArraysUL - 1);
      if(!settings.pAltAlns)
      { /*If: I am not printing full alt alignments*/
         pAltAlnScores(
            alnMtrxST,
            settings.minScoreL,
            altAlnFILE
         ); /*Print out socres for all alignments*/

         if(settings.justScoresBl)
         { /*If: I am printing out scores only*/
            if(outFILE != stdout) fclose(outFILE);
            if(altAlnFILE != stdout) fclose(altAlnFILE);

            freeSeqSTStack(&refST);
            freeSeqSTStack(&queryST);
            freeAlnMatrix(alnMtrxST);
   
            exit(0);
         } /*If: I am printing out scores only*/
      } /*If: I am not printing full alt alignments*/

      else
      { /*Else: printing out an alignment for each alt*/
         for(
            iterUL = 0;
            iterUL < alnMtrxST->lenArraysUL;
            ++iterUL
         ){ /*Loop: Pint out kept alignments*/

            indexToCoord(
               alnMtrxST->lenRefUL,
               alnMtrxST->startIndexAryUL[iterUL],
               refST.offsetUL,
               queryST.offsetUL
            );
      
            indexToCoord(
               alnMtrxST->lenRefUL,
               alnMtrxST->endIndexAryUL[iterUL],
               refST.endAlnUL,
               queryST.endAlnUL
            );

            alnST = Hirschberg(&refST,&queryST,&settings);
            if(alnST == 0) goto memWaterAltErr;

            errUC =             
               printAln(
                  altAlnFILE,
                  altAlnFileStr,
                  &refST,
                  &queryST,
                  alnST,
                  alnMtrxST->scoreAryL[iterUL],
                  &settings,
                  scoreMtrxFileStr/*Name: scoring matrix*/
            ); /*Print out the alternative alignment*/

            freeAlnST(alnST);
            alnST = 0;

            if(errUC)
            { /*If: I falied to print out the alignment*/
               memWaterAltErr:

               if(altAlnFILE != stdout)
                  fclose(altAlnFILE);

               if(outFILE != stdout) fclose(outFILE);

               freeSeqSTStack(&refST);
               freeSeqSTStack(&queryST);
               freeAlnMatrix(alnMtrxST);
   
               fprintf(
                  stderr,
                  "Error while printing alternative"
               );
               fprintf(stderr,  "alignments\n");

               exit(-1);
            } /*If: I falied to print out the alignment*/
         } /*Loop: Pint out kept alignments*/
      } /*Else: printing out an alignment for each alt*/

      indexToCoord(
         alnMtrxST->lenRefUL,
         alnMtrxST->bestStartIndexUL,
         refST.offsetUL,
         queryST.offsetUL
      );

      indexToCoord(
         alnMtrxST->lenRefUL,
         alnMtrxST->bestEndIndexUL,
         refST.endAlnUL,
         queryST.endAlnUL
      );

      alnST = Hirschberg(&refST, &queryST, &settings);

      bestScoreL = alnMtrxST->bestScoreL;
      freeAlnMatrix(alnMtrxST); /*No longer need*/
      alnMtrxST = 0;

      goto printAlignment;
   } /*Else if; memory waterman alignment + scan*/

   /*****************************************************\
   * Main Sec-05 Sub-05:
   *  - Check if doing an memory efficent Waterman
   *    alignment that prints out a single alignment
   *    (uses a Hirschberg to get the alignment)
   \*****************************************************/

   else if(settings.memWaterBl)
   { /*Else I am just finding the best alignment*/
     if(settings.noGapBl)
        alnMtrxST =
           memWaterNoGap(&queryST,&refST,&settings);
     else alnMtrxST = memWater(&queryST,&refST,&settings);

      if(alnMtrxST == 0)
      { /*If: the aligment falied*/
         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);
       
         fprintf(
            stderr,
            "Ran out of memory for mem-water alignment\n"
         );
       
         if(outFILE != stdout) fclose(outFILE);
         if(altAlnFILE != stdout) fclose(altAlnFILE);
   
         freeSeqSTStack(&refST);
         freeSeqSTStack(&queryST);
   
         exit(-1);
      } /*If: the aligment falied*/

      indexToCoord(
         alnMtrxST->lenRefUL,
         alnMtrxST->bestStartIndexUL,
         refST.offsetUL,
         queryST.offsetUL
      );

      indexToCoord(
         alnMtrxST->lenRefUL,
         alnMtrxST->bestEndIndexUL,
         refST.endAlnUL,
         queryST.endAlnUL
      );

      if(settings.justScoresBl)
      { /*If I am just printing out coordinates*/
         fprintf(
            outFILE,
            "%li\t%lu\t%lu\t%lu\t%lu\n",
            alnMtrxST->bestScoreL,
            refST.offsetUL,
            refST.endAlnUL,
            queryST.offsetUL,
            queryST.endAlnUL
         );

         goto noAlnOutFree;
      } /*If I am just printing out coordinates*/

      alnST = Hirschberg(&refST, &queryST, &settings);

      bestScoreL = alnMtrxST->bestScoreL;
      freeAlnMatrix(alnMtrxST); /*No longer need*/
      alnMtrxST = 0;

      goto printAlignment;
   } /*Else I am just finding the best alignment*/

   /******************************************************\
   * Main Sec-05 Sub-06:
   *  - Check if doing an Hirschberg alignment
   \******************************************************/

   else if(settings.useHirschBl != 0)
   { /*Else if doing an Hirschberg alignment*/
     if(settings.noGapBl)
       alnST = HirschbergNoGap(&refST,&queryST,&settings);
     else
        alnST = Hirschberg(&refST, &queryST, &settings);

     goto printAlignment;
     /* The Hirschberg returns an alignment structure,
     `  instead of an directional matrix.
     */
   } /*Else if doing an Hirschberg alignment*/

   /******************************************************\
   * Main Sec-05 Sub-07:
   *  - Let user know this was an invalid alignment
   \******************************************************/
   else
   { /*If no aignment was requested*/
      freeSeqSTStack(&refST);
      freeSeqSTStack(&queryST);

      printHelpMesg(stderr, 1); /*short help*/
      fprintf(
         stderr,
         "Invalid or no aligment method input\n"
      ); /*Let user know the error*/

      exit(-1);
   } /*If no aignment was requested*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-06:
   ^  - Print out alignments and clean up
   ^  o main sec-06 sub-01:
   ^    - Check if have an alignment structure to print
   ^      out
   ^  o main sec-06 sub-02:
   ^    - Print out the alignments
   ^  o main sec-06 sub-03:
   ^    - Clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   *  Main Sec-06 Sub-01:
   *   - Check if have an alignment structure to print
   *     out
   \*****************************************************/

   printAlignment:

   if(alnST == 0)
   { /*If: I falied to make an alignment array*/
      if(outFILE != stdout) fclose(outFILE);
      if(altAlnFILE != stdout) fclose(altAlnFILE);

      freeSeqSTStack(&refST);
      freeSeqSTStack(&queryST);

       fprintf(
         stderr,
         "Memory error when making alignment array\n"
       );

       exit(-1);
   } /*If: I falied to make an alignment array*/

   /*****************************************************\
   *  Main Sec-06 Sub-02:
   *   - Print out the alignments
   \*****************************************************/

   lookupIndexToSeq(refST.seqCStr);
   lookupIndexToSeq(queryST.seqCStr);

   if(
      printAln(
         outFILE,
         outFileStr,
         &refST,
         &queryST,
         alnST,
         bestScoreL,
         &settings,
         scoreMtrxFileStr /*File name for scoring matrix*/
      )
   ){ /*If could not print out the alignmnet*/
      fprintf(stderr, "Failed to print alignment\n");

      if(outFILE != stdout) fclose(outFILE);
      if(altAlnFILE != stdout) fclose(altAlnFILE);

      freeAlnST(alnST); /*NEED TO SET UP*/
      freeSeqSTStack(&refST);
      freeSeqSTStack(&queryST);

      exit(1);
   } /*If could not print out the alignmnet*/

   /*****************************************************\
   *  Main Sec-06 Sub-03:
   *   - Clean up and exit
   \*****************************************************/

   freeAlnST(alnST);

   noAlnOutFree: /*When memWater just printing positions*/

   if(outFILE != stdout) fclose(outFILE);
   if(altAlnFILE != stdout) fclose(altAlnFILE);

   freeSeqSTStack(&refST);
   freeSeqSTStack(&queryST);

   exit(0);
} /*main*/

/*-------------------------------------------------------\
| Fun-01: checkInput
|  - Gets user input and lets program now if something
|    was incorrect.
| Input:
|  - lenArgsInt:
|    o Number of arguments the user input
|  - argsCStr:
|    o Array of c-strings with the users input arguments
|  - refFileCStr:
|    o Pointer to c-string to hold the name of the fasta
|      file with the reference sequence
|  - qryFileCStr:
|    o Pointer to c-string to hold the name of the fasta
|      file with the query sequence
|  - outFileStr:
|    o Pointer to c-string to hold the name of the file to
|      output the main alignment to
|  - altAlnFileStr:
|    o Pointer to c-string to hold the name of the file to
|      output the alternative alignments to
|  - refStartAlnUL
|    o Pointer to unsigned long to hold the first position
|      to align on the reference sequence
|  - refEndAlnUL
|    o Pointer to unsigned long to hold the last position
|      to align on the reference sequence
|  - qryStartAlnUL
|    o Pointer to unsigned long to hold the first position
|      to align on the query sequence
|  - qryEndAlnUL
|    o Pointer to unsigned long to hold the last position
|      to align on the query sequence
|  - scoreMtrxFileStr:
|    o Pointer to hold the file name of the input score
|      matrix (socre matrix is read into settings by
|      checkInput)
|  - matchMtrxFileStr:
|    o Pointer to hold the file name of the input match
|      matrix (the match matrix is read into settings by
|      checkInput)
|  - settings:
|    o Pointer to alnSet structure to hold the users
|      input settings for the alignment
| Output:
|  - Modifies:
|    o Each input variable to hold the specified user
|      input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to "-score-matrix" for invalid scoring
|       matrix input
|  - Prints to stdout when the scoring file is invalid
\-------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,     /*Number arguments user input*/
    char *argsCStr[],    /*Array with user arguments*/
    char **refFileCStr,  /*file name of reference file*/
    char **queryFileCStr,/*File name of the query file*/
    char **outFileStr,     /*Name of the output file*/
    char **altAlnFileStr, /*file: alternative alignments*/
    ulong *refStartAlnUL, /*Start of reference alignment*/
    ulong *refEndAlnUL,   /*End of reference alignment*/
    ulong *qryStartAlnUL, /*Start of query alignment*/
    ulong *qryEndAlnUL,   /*End of query alignment*/
    char **scoreMtrxFileStr,/*Holds scoring matrix file*/
    char **matchMtrxFileStr,/*Holds matching matrix file*/
    struct alnSet *settings /*Aligment settings*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '  o fun-01 sec-01:
   '    -  Variable declerations
   '  o fun-01 sec-02:
   '    - Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-1:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    char *singleArgCStr = 0;
    unsigned long fileErrUL = 0;
    FILE *inFILE = 0; /*For loading the scoring matrix*/

    int iArg = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-01 Sec-02:
    ^  - Look through user input
    ^  o fun-01 sec-02 sub-01:
    ^    - Loop though user input and get parm/arg
    ^  o fun-01 sec-02 sub-02:
    ^    - Required variables (reference/query fasta)
    ^  o fun-01 sec-02 sub-03:
    ^    - Choosing the alignment algorithm
    ^  o fun-01 sec-02 sub-04:
    ^    - Main file format options
    ^  o fun-01 sec-02 sub-05:
    ^    - alignment file output format options
    ^  o fun-01 sec-02 sub-06:
    ^    - Query reference scan options
    ^  o fun-01 sec-02 sub-07:
    ^    - General alignment options
    ^  o fun-01 sec-02 sub-08:
    ^    - Read in the scoring matrix
    ^  o fun-01 sec-02 sub-09:
    ^    - Read in the match matrix
    ^  o fun-01 sec-02 sub-10:
    ^    - Report errors or when finshed exit
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /****************************************************\
    * Fun-01 Sec-02 sub-01:
    *  - Loop though user input and get parm/arg
    \****************************************************/

    for(iArg = 1; iArg < *lenArgsInt; ++iArg)
    { /*loop through all user input arguments*/
        /*The 0 index holds the program name*/
        singleArgCStr = *(argsCStr +iArg + 1);/*arg*/
        tmpCStr = *(argsCStr + iArg);       /*Paramter*/

       /*************************************************\
       * Fun-01 Sec-02 sub-02:
       *  - Required variables (reference/query fasta)
       \*************************************************/

        if(strcmp(tmpCStr, "-ref") == 0)
        { /*If: the user provied a reference sequence*/
            *refFileCStr = singleArgCStr;
            ++iArg;
        } /*If: the user provied a reference sequence*/

        else if(strcmp(tmpCStr, "-query") == 0)
        { /*If: the user provied a query sequence*/
            *queryFileCStr = singleArgCStr;
            ++iArg;
        } /*If: the user provied a query sequence*/

       /*************************************************\
       * Fun-01 Sec-02 sub-03:
       *  - Choosing the alignment algorithm
       \*************************************************/

       else if(strcmp(tmpCStr, "-use-needle") == 0)
       { /*Else if using a needleman alignment*/
           settings->useNeedleBl = 1;
           settings->useWaterBl = 0;
           settings->useHirschBl = 0;
           settings->memWaterBl = 0;
       } /*Else if using a needleman alignment*/

       else if(strcmp(tmpCStr, "-use-water") == 0)
       { /*Else if doing a waterman smith alignment*/
           settings->useNeedleBl = 0;
          settings->useWaterBl = 1;
           settings->useHirschBl = 0;
           settings->memWaterBl = 0;
       } /*Else if doing a waterman smith alignment*/

       else if(strcmp(tmpCStr, "-use-hirschberg") == 0)
       { /*Else if doing a Hirshberg alignment*/
           settings->useNeedleBl = 0;
           settings->useWaterBl = 0;
           settings->useHirschBl = 1;
           settings->memWaterBl = 0;
       } /*Else if doing a Hirshberg alignment*/

       else if(strcmp(tmpCStr, "-use-mem-water") == 0)
       { /*Else if I am doing a memory effecient water*/
           settings->useNeedleBl = 0;
           settings->useWaterBl = 0;
           settings->useHirschBl = 0;
           settings->memWaterBl = 1;
       } /*Else if I am doing a memory effecient water*/

       /*************************************************\
       * Fun-01 Sec-02 sub-04:
       *  - Main file format options
       \*************************************************/

        /*Output file*/
        else if(strcmp(tmpCStr, "-out") == 0)
        { /*Else If: Have an output file*/
            *outFileStr = singleArgCStr;
            ++iArg;
        } /*Else If: Have an output file*/

        /*Alternative alignment output file*/
        else if(strcmp(tmpCStr, "-alt-out") == 0)
        { /*Else If: Alt alignments in separate file*/
            *altAlnFileStr = singleArgCStr;
            ++iArg;
        } /*Else If: Alt alignments in separate file*/

        /*Only print out the aligned regions*/
        else if(strcmp(tmpCStr, "-print-aligned") == 0)
           settings->pFullAlnBl = 0;

        /*Printing out unligned (softmasked) regions*/
        else if(strcmp(tmpCStr, "-print-unaligned") == 0)
           settings->pFullAlnBl = 1;

        /*Printing out base positions for each line*/
        else if(strcmp(tmpCStr, "-print-positions") == 0)
           settings->pBasePosBl = 1;

        /*Not printing out base positions for each line*/
        else if(strcmp(tmpCStr, "-no-positions") == 0)
           settings->pBasePosBl = 0;

        else if(strcmp(tmpCStr, "-line-wrap") == 0)
        { /*Else If: the user provided the line wrap*/
           base10StrToUI(
              singleArgCStr,
              settings->lineWrapUS
           );

           ++iArg;
        } /*Else If: the user provided the line wrap*/

       /*************************************************\
       * Fun-01 Sec-02 sub-05:
       *  - alignment file output format options
       \*************************************************/

        /*using the expanded cigar format*/
        else if(strcmp(tmpCStr, "-format-expand-cig")==0)
           settings->formatFlag = defExpandCig;

        /*using an EMBOSS like format*/
        else if(strcmp(tmpCStr, "-format-emboss") == 0)
           settings->formatFlag = defEMBOSS;

        /*using clustal format*/
        else if(strcmp(tmpCStr, "-format-clustal") == 0)
           settings->formatFlag = defClustal;

        /*Using fasta format*/
        else if(strcmp(tmpCStr, "-format-fasta") == 0)
           settings->formatFlag = defFasta;

       /*************************************************\
       * Fun-01 Sec-02 sub-06:
       *  - Query reference scan options
       \*************************************************/

        else if(
          strcmp(tmpCStr, "-query-ref-scan") == 0
        ) { /*Else If: doing alternative alignments*/
          settings->refQueryScanBl = 1;
          settings->useNeedleBl = 0;
          settings->useHirschBl = 0;

          if(!(settings->useWaterBl))
             settings->memWaterBl = 1;
        } /*Else If: doing alternative alignments*/

        /*Do not print out the main alignment*/
        else if(strcmp(tmpCStr, "-only-scores") == 0)
           settings->justScoresBl = 1;

        /*Print out the main aligmnent with and scores*/
        else if(strcmp(tmpCStr, "-scores-and-aln") == 0)
           settings->justScoresBl = 0;

        else if(strcmp(tmpCStr, "-min-score") == 0)
        { /*Else If: have a min score to keep alignment*/
           base10StrToSL(
              singleArgCStr,
              settings->minScoreL
           );

           ++iArg;
        } /*Else If: have a min score to keep alignment*/

        /*Check if filtering by query overlaps*/
        else if(strcmp(tmpCStr, "-filt-query") == 0)
           settings->filtQryBl = 1;

        else if(strcmp(tmpCStr, "-no-filt-query") == 0)
           settings->filtQryBl = 0;

        /*Check if filtering by reference overlaps*/
        else if(strcmp(tmpCStr, "-filt-ref") == 0)
           settings->filtRefBl = 1;

        else if(strcmp(tmpCStr, "-no-filt-ref") == 0)
           settings->filtRefBl = 0;

        /*if filtering by reference and query overlaps*/
        else if(strcmp(tmpCStr, "-filt-ref-query")==0)
           settings->filtQryRefBl = 1;

        else if(strcmp(tmpCStr, "-no-filt-ref-query")==0)
           settings->filtQryRefBl = 0;

        /*Print out full alternative alignments*/
        else if(strcmp(tmpCStr, "-p-alt-alignments") == 0)
           settings->pAltAlns = 1;

       /*************************************************\
       * Fun-01 Sec-02 sub-07:
       *  - General alignment options
       \*************************************************/

       else if(strcmp(tmpCStr, "-no-gapextend") == 0)
          settings->noGapBl = 1;
       else if(strcmp(tmpCStr, "-use-gapextend") == 0)
          settings->noGapBl = 0;

       else if(strcmp(tmpCStr, "-gapopen") == 0)
       { /*Else If: getting the gap opening penalty*/
         base10StrToSC(singleArgCStr, settings->gapOpenC);
         ++iArg;
       } /*Else If: getting the gap opening penalty*/

       else if(strcmp(tmpCStr, "-gapextend") == 0)
       { /*Else If: getting the gap extension penalty*/
         settings->noGapBl = 0;
         base10StrToSC(
            singleArgCStr,
            settings->gapExtendC
         );
         ++iArg;
       } /*Else If: getting the gap extension penalty*/

       else if(strcmp(tmpCStr, "-two-bit") == 0)
          settings->twoBitBl = 1;

       else if(strcmp(tmpCStr, "-no-two-bit") == 0)
          settings->twoBitBl = 0;

       else if(strcmp(tmpCStr, "-ref-start") == 0)
       { /*Else If: start of reference alignment input*/
          base10StrToUL(singleArgCStr, *refStartAlnUL);
          --(*refStartAlnUL); /*Convert to index 0*/
          ++iArg;
       } /*Else If: start of reference alignment input*/

       else if(strcmp(tmpCStr, "-ref-end") == 0)
       { /*Else If: end of reference alignment input*/
          base10StrToUL(singleArgCStr, *refEndAlnUL);
          --(*refEndAlnUL); /*Convert to index 0*/
          ++iArg;
       } /*Else If: end of reference alignment input*/

       else if(strcmp(tmpCStr, "-query-start") == 0)
       { /*Else If: start of query alignment input*/
          base10StrToUL(singleArgCStr, *qryStartAlnUL);
          --(*qryStartAlnUL); /*Convert to index 0*/
          ++iArg;
       } /*Else If: start of query alignment input*/

       else if(strcmp(tmpCStr, "-query-end") == 0)
       { /*Else If: end of query alignment input*/
          base10StrToUL(singleArgCStr, *qryEndAlnUL);
          --(*qryEndAlnUL); /*Convert to index 0*/
          ++iArg;
       } /*Else If: end of query alignment input*/

        /*Only do these checks when the user has not
        ` called a direction flag during compile time
        */
        #if defined SNPINSDEL
        #elif defined SNPDELINS
        #elif defined INSSNPDEL 
        #elif defined INSDELSNP
        #elif defined DELSNPINS
        #elif defined DELINSSNP
        #else
           else if(strcmp(tmpCStr, "-snp-ins-del") == 0)
               settings->bestDirC = defSnpInsDel;

           else if(strcmp(tmpCStr, "-snp-del-ins") == 0)
               settings->bestDirC = defSnpDelIns;

           else if(strcmp(tmpCStr, "-ins-snp-del") == 0)
               settings->bestDirC = defInsSnpDel;

           else if(strcmp(tmpCStr, "-del-snp-ins") == 0)
               settings->bestDirC = defDelSnpIns;

           else if(strcmp(tmpCStr, "-ins-del-snp") == 0)
               settings->bestDirC = defInsDelSnp;

           else if(strcmp(tmpCStr, "-del-ins-snp") == 0)
               settings->bestDirC = defDelInsSnp;
        #endif

       /*************************************************\
       * Fun-01 Sec-02 sub-08:
       *  - Read in the scoring matrix
       \*************************************************/

        else if(strcmp(tmpCStr, "-score-matrix") == 0)
        { /*else if the user supplied a scoring matrix*/
            *scoreMtrxFileStr = singleArgCStr;
            inFILE = fopen(singleArgCStr, "r");

            if(inFILE == 0)
            { /*If I could not open the scoring file*/
                return tmpCStr; /*So user knows invalid*/

                fprintf(stderr,
                  "-score-matrix %s is not an file\n",
                  singleArgCStr
                ); /*Print out the problem*/
            } /*If I could not open the scoring file*/

            fileErrUL = readInScoreFile(settings, inFILE);

            if(fileErrUL != 0)
            { /*If the scoring file had an errors*/
              fprintf(
               stderr,
               "Invalid line (%lu) in -score-matrix %s\n",
               fileErrUL,
               singleArgCStr
              ); /*Print out the problem*/

                return tmpCStr;    /*invalid file*/
            } /*If the scoring file had an error*/

            ++iArg;
        } /*else if the user supplied a scoring matrix*/

       /*************************************************\
       * Fun-01 Sec-02 sub-09:
       *  - Read in the match matrix
       \*************************************************/

        else if(strcmp(tmpCStr, "-match-matrix") == 0)
        { /*else if the user supplied a scoring matrix*/
            *matchMtrxFileStr = singleArgCStr;
            inFILE = fopen(singleArgCStr, "r");

            if(inFILE == 0)
            { /*If I could not open the scoring file*/
                return tmpCStr; /*So user knows invalid*/

                fprintf(stderr,
                  "-match-matrix %s is not an file\n",
                  singleArgCStr
                ); /*Print out the problem*/
            } /*If I could not open the scoring file*/

            fileErrUL = readInMatchFile(settings, inFILE);

            if(fileErrUL != 0)
            { /*If the scoring file had an errors*/
              fprintf(
               stderr,
               "Invalid line (%lu) in -match-matrix %s\n",
               fileErrUL,
               singleArgCStr
              ); /*Print out the problem*/

                return tmpCStr;    /*invalid file*/
            } /*If the scoring file had an error*/

            ++iArg;
        } /*else if the user supplied a scoring matrix*/

       /*************************************************\
       * Fun-01 Sec-02 sub-10:
       *  - Report errors or when finshed exit
       \*************************************************/
            
        else return tmpCStr; /*Invalid parameter*/
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/

/*-------------------------------------------------------\
| Fun-02: printHelpMesg
|  - Prints the help message for alnSeq
| Input:
|  - outFILE:
|    o File to print the help message to
|  - breifBl:
|    o 1: Print the short help message
|    o 0: Print the entire help message
| Output:
|  - Prints:
|    o help message to outFILE.
\-------------------------------------------------------*/
void printHelpMesg(
   FILE *outFILE,
   char breifBl   /*Print a shorter help message*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: printHelpMesg
   '  - Prints out the help message to outFILE
   '  o fun-02 Sec-01:
   '    - Usage block
   '  o fun-02 Sec-02:
   '    - Input block
   '  o fun-02 sec-03:
   '    - Output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-01:
   ^  - Usage block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      outFILE,
      "alnSeq -query query.fasta -ref ref.fasta"
   );
   fprintf(outFILE, " [options...]\n");

   fprintf(
      outFILE,
      "  - Does a pairwise alignment on two input"
   );
   fprintf(outFILE, "sequences \n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Input block
   ^  o fun-02 sec-02 sub-01:
   ^     - General input/output block (includes line wrap)
   ^  o fun-02 sec-02 sub-02:
   ^     - Alignment algorithim selection block
   ^  o fun-02 sec-02 sub-03:
   ^     - Alignment paramaters block
   ^  o fun-02 sec-02 sub-04:
   ^    - File output settings (non-format)
   ^  o fun-02 sec-02 sub-05:
   ^     - File output format block
   ^  o fun-02 sec-02 sub-06:
   ^     - Waterman specific paramters block
   ^  o fun-02 sec-02 sub-07:
   ^     - Selecting alignment direction block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-01:
   *  - General input/output block (includes line wrap)
   \******************************************************/

   fprintf(outFILE, "Input:\n");
   fprintf(outFILE, "  General:\n");

   fprintf(outFILE, "    -query: [Required]\n");
   fprintf(
      outFILE,
      "      o Fasta file with query sequence.\n"
   );

   fprintf(outFILE, "    -ref: [Required]\n");
   fprintf(
      outFILE,
      "      o Fasta file with reference sequence.\n"
   );

   fprintf(outFILE, "    -out: [stdout]\n");
   fprintf(
      outFILE,
      "      o File to output alignment to"
   );
   fprintf(outFILE, " (use \"-\" for stdout).\n");



   /******************************************************\
   * Fun-02 Sec-02 Sub-02:
   *  - Alignment algorithim selection block
   \******************************************************/

   fprintf(outFILE, "  Algorithims:\n");
   if(defUseNeedle)
      fprintf(outFILE, "    -use-needle: [Yes]\n");
   else fprintf(outFILE, "    -use-needle: [No]\n");

   fprintf(
      outFILE,
      "      o Do a Needleman Wunsch alignment.\n" 
   );

   if(defUseWater)
      fprintf(outFILE, "    -use-water: [Yes]\n");
   else fprintf(outFILE, "    -use-water: [No]\n");

   fprintf(
      outFILE,
      "      o Do a Waterman Smith alignment.\n"
   );

   if(defUseHirsch)
       fprintf(outFILE, "    -use-hirschberg: [Yes]\n");
   else fprintf(outFILE, "    -use-hirschberg: [No]\n");

   fprintf(
      outFILE,
      "      o Do a Hirschberg alignment.\n"
   );

   if(defUseMemWater)
       fprintf(outFILE, "    -use-mem-water: [Yes]\n");
   else fprintf(outFILE, "    -use-mem-water: [No]\n");

   fprintf(
      outFILE,
      "      o Do a memory efficent Waterman alignment.\n"
   );

   fprintf(
     outFILE,
     "      o Uses a Waterman to find the starting and"
   );
   fprintf(outFILE, " ending\n        positions.");
   fprintf(
      outFILE,
      " A Hirschberg is then used to find\n"
   );
   fprintf(outFILE, "        the alignment.\n");

   /******************************************************\
   * Fun-02 Sec-02 Sub-03:
   *  - Alignment paramaters block
   \******************************************************/

   fprintf(outFILE, "  Alignment settings:\n");

   if(defUseTwoBit)
      fprintf(outFILE, "    -two-bit: [Yes]\n");
   else
      fprintf(outFILE, "    -two-bit: [No]\n");

   fprintf(
      outFILE,
      "      o Use two bit arrays in an alignment."
   );
   fprintf(outFILE, " Takes less\n");
   fprintf(
      outFILE,
      "        memory (~1/4), but doubles the time.\n"
   );

   fprintf(
      outFILE,
      "      o This is ignored for all memWater and"
   );
   fprintf(outFILE, " Hirschberg\n");
   fprintf(
      outFILE,
      "        alignments.\n"
   );
   fprintf(
      outFILE,
      "      o Disable with -no-two-bit\n"
   );


   fprintf(outFILE, "    -gapopen: [%i]\n", defGapOpen);
   fprintf(
      outFILE,
      "      o Cost of starting an indel (as integer). A "
   );
   fprintf(
      outFILE,
      "negative\n        value is a penalty.\n"
   );

   fprintf(outFILE,"    -gapextend: [%i]\n",defGapExtend);
   fprintf(
     outFILE,
     "      o Cost of extending an indel one base. A"
   );
   fprintf(
      outFILE,
      " number < 0\n        is a penalty.\n"
   );
   fprintf(
     outFILE,
     "      o Disable: -no-gapextend;"
   );
   fprintf(
      outFILE,
      " Enable: -use-gapextend\n"
   );

   fprintf(
      outFILE,
      "    -score-matrix: [%s]\n",
      defMatrixNameStr
   );
   fprintf(
      outFILE,
      "      o File with scoring matrix to use. It should"
   );
   fprintf(
     outFILE,
     " have\n        one line for each score to change.\n"
   );
   fprintf(
      outFILE,
      "      o (see scoring-matrix.txt for an example).\n"
   );

   fprintf(outFILE, "    -match-matrix: [DNA]\n");
   fprintf(
      outFILE,
      "      o File with match matrix to use. It has one"
   );
   fprintf(
     outFILE,
     " line\n        for each possible base pair.\n"
   );
   fprintf(
      outFILE,
      "      o (see match-matrix.txt for an example).\n"
   );
 
   if(breifBl) goto helpOutputBlock;

   fprintf(outFILE, "    -ref-start: 1\n");
   fprintf(
      outFILE,
      "      o First base in reference align (index 1)\n"
   );

   fprintf(outFILE, "    -ref-end: length(reference)\n");
   fprintf(
      outFILE,
      "      o Last base in reference align (index 1)\n"
   );

   fprintf(outFILE, "    -query-start: 1\n");
   fprintf(
      outFILE,
      "      o First base in query align (index 1)\n"
   );

   fprintf(outFILE, "    -query-end: length(query)\n");
   fprintf(
      outFILE,
      "      o Last base in query align (index 1)\n"
   );

   /*****************************************************\
   * Fun-02 Sec-02 Sub-04:
   *  - Query reference scan settings run settings
   \*****************************************************/

   fprintf(outFILE, "  Query reference scan:\n");
   fprintf(
      outFILE,
      "    o A query reference scan saves a score for"
   );
   fprintf(outFILE, "each\n");
   fprintf(
      outFILE,
      "      reference and query base if possible.\n"
   );

   if(defQueryRefScan)
       fprintf(outFILE,"    -query-ref-scan: [Yes]\n");
   else fprintf(outFILE,"    -query-ref-scan: [No]\n");

   fprintf(
      outFILE,
      "      o prints out the best score for each"
   );
   fprintf(
      outFILE,
      " reference and\n        query base.\n"
   );
   fprintf(
      outFILE,
      "      o Options are -use-mem-water [Default] or.\n"
   );
   fprintf(outFILE, "        -use-water.\n");

   fprintf(outFILE, "    -alt-out: [file from -out]\n");
   fprintf(
     outFILE,
     "      o File to output alternative alignments to.\n"
   );
   fprintf(outFILE, "      o use \"-\" for stdout.\n");


   fprintf(outFILE,"    -min-score: [%i]\n", defMinScore);
   fprintf(
      outFILE,
      "      o Minimum score needed to keep an\n"
   );
   fprintf(outFILE, "        alternative alignment.\n");


   if(defFilterByQueryRef)
       fprintf(outFILE,"    -filt-ref-query: [Yes]\n");
   else fprintf(outFILE,"    -filt-ref-query: [No]\n");

   fprintf(
      outFILE,
      "      o Remove alternative alignments that overlap"
   );

   fprintf(
      outFILE,
      "\n        in both the reference and query.\n"
   );
   fprintf(
      outFILE,
      "      o Disable with -no-filt-ref-query).\n"
   );

   if(defFilterByRef)
       fprintf(outFILE,"    -filt-ref: [Yes]\n");
   else fprintf(outFILE,"    -filt-ref: [No]\n");

   fprintf(
      outFILE,
      "      o Remove alternative alignments that overlap"
   );

   fprintf(
      outFILE,
      "\n        in the reference, (ignores overlaps in"
   );
   fprintf(outFILE, " query).\n");

   fprintf(
      outFILE,
      "      o Disable with -no-filt-ref).\n"
   );

   if(defFilterByQuery)
       fprintf(outFILE,"    -filt-query: [Yes]\n");
   else fprintf(outFILE,"    -filt-query: [No]\n");

   fprintf(
      outFILE,
      "      o Remove alternative alignments that overlap"
   );

   fprintf(
      outFILE,
      "\n        in the query, (ignores overlaps in"
   );
   fprintf(outFILE, " reference).\n");
   fprintf(
      outFILE,
      "      o Disable with -no-filt-query).\n"
   );

   /*****************************************************\
   * Fun-02 Sec-02 Sub-05:
   *  - Query reference scan output settings
   \*****************************************************/

   fprintf(outFILE, "  Query reference scan output:\n");

   fprintf(outFILE, "    -alt-out: [-out]\n");
   fprintf(
      outFILE,
      "      o Print alternative alignments to this file."
   );
   fprintf(outFILE, "\n");

   if(defJustScoresBl)
      fprintf(outFILE,"    -only-scores: [Yes]\n");
   else
      fprintf(outFILE,"    -only-scores: [No]\n");

   fprintf(
      outFILE,
      "      o Only print out the score, starting"
   );
   fprintf(outFILE, " coordinates\n");
   fprintf(
      outFILE,
      "        and ending coordinaates of an alignment\n"
   );
   fprintf(
      outFILE,
      "      o Only for the memory effiecent Watermans\n"
   );
   fprintf(
      outFILE,
      "      o Disable with -scores-and-aln\n"
   );

   if(defPAltAln) 
      fprintf(outFILE,"    -p-alt-alignments: [Yes]\n");
   else
      fprintf(outFILE,"    -p-alt-alignments: [No]\n");

   fprintf(
     outFILE,
     "      o Print out an alignment for each alternate\n"
   );
   fprintf( 
      outFILE,
      "        alignment in a query reference scan\n"
   );

   /*****************************************************\
   * Fun-02 Sec-02 Sub-06:
   *  - General output formatting
   \*****************************************************/

   fprintf(outFILE, "  General output formating:\n");
   if(defPAln)
      fprintf(outFILE,"    -print-aligned: [Yes]\n");
   else fprintf(outFILE,"    -print-aligned: [No]\n");

   fprintf(
      outFILE,
      "      o Only print the aligned regions of each"
   );
   fprintf(outFILE, " sequence.\n");

   if(!defPAln)
      fprintf(outFILE,"    -print-unaligned: [Yes]\n");
   else fprintf(outFILE,"    -print-unaligned: [No]\n");

   fprintf(
      outFILE,
      "      o Print out the entire reference and query"
   );
   fprintf(outFILE, " sequence.\n");

   if(defPPos)
      fprintf(outFILE, "    -print-positions: [Yes]\n");
   else fprintf(outFILE, "    -print-positions: [No]\n");

   fprintf(
      outFILE,
      "      o Print the starting and ending position for"
   );
   fprintf(outFILE, " each\n");
   fprintf(outFILE, "     line in the alignment.\n");
   fprintf(
      outFILE,
     "      o Clustal format only prints ending position"
   );
   fprintf(outFILE, ".\n");

   if(!defPPos)
      fprintf(outFILE, "    -no-positions: [Yes]\n");
   else fprintf(outFILE, "    -no-positions: [No]\n");

   fprintf(
      outFILE,
      "      o Turns off -print-positions\n"
   );
   fprintf(
      outFILE,
      "      o EMBOSS format always prints positions.\n"
   );


   fprintf(outFILE,"    -line-wrap: [%i]\n", defLineWrap);
   fprintf(
    outFILE,
    "        o Number characters per line in output file"
   );
   fprintf(outFILE, ".\n");
   fprintf(outFILE, "      o Use 0 for no line wrap.\n");
   fprintf(
      outFILE,
      "      o Minimum line wrap is 10 for fasta, 32 for"
   );
   fprintf(outFILE, "\n");
   fprintf(
      outFILE,
      "        clustal, and 42 for expanded cigar/EMBOSS"
   );
   fprintf(outFILE, "\n");

   /*****************************************************\
   * Fun-02 Sec-02 Sub-07:
   *  - File output format block
   \*****************************************************/

   fprintf(outFILE, "  File formats:\n");
   if(defFormat == defExpandCig)
      fprintf(outFILE, "    -format-expand-cig: [Yes]\n");
   else
      fprintf(outFILE, "    -format-expand-cig: [No]\n");

   fprintf(
      outFILE,
      "      o Prints the reference sequence, then query"
   );
   fprintf(outFILE, "\n        sequence,");
   fprintf(outFILE, " and then the eqx line.\n");
   fprintf(outFILE, "      o eqx line format:\n");
   fprintf(outFILE, "        - I = Insertion\n");
   fprintf(outFILE, "        - D = Deletion\n");
   fprintf(outFILE, "        - X = mismatch\n");
   fprintf(outFILE, "        - = = match\n");
   fprintf(outFILE, "       - S = soft mask\n");
   fprintf(outFILE, "      o reference/query lines\n");
   fprintf(
      outFILE,
      "        - Without -print-positions\n"
   );

   fprintf(outFILE, "          - Ref: sequence\n");
   fprintf(outFILE, "          - Qry: sequence\n");
   fprintf(
      outFILE,
      "        - With -print-positions\n"
   );
   fprintf(
      outFILE,
      "          - Ref: start-base sequence end-base\n"
   );
   fprintf(
      outFILE,
      "          - Qry: start-base sequence end-base\n"
   );


   if(defFormat == defEMBOSS)
      fprintf(outFILE, "    -format-emboss: [Yes]\n");
   else fprintf(outFILE, "    -format-emboss: [No]\n");

   fprintf(
      outFILE,
      "      o Prints the reference sequence, then"
   );
   fprintf(outFILE," eqx line, and\n");
   fprintf(
       outFILE,
       "        then the query sequence.\n"
   );
   fprintf(outFILE, "      o eqx line format:\n");
   fprintf(outFILE, "        - | = Match\n");
   fprintf(outFILE, "        - space = SNP/gap\n");
   fprintf(outFILE, "      o reference/query lines\n");
   fprintf(
      outFILE,
      "        - ID: start-base sequence end-base\n"
   );

   if(defFormat == defClustal)
      fprintf(outFILE, "    -format-clustal: [Yes]\n");
   else fprintf(outFILE, "    -format-clustal: [No]\n");

   fprintf(
      outFILE,
      "      o Prints the reference sequence, then query"
   );
   fprintf(outFILE, "\n        sequence,");
   fprintf(outFILE, " and then the eqx line.\n");
   fprintf(outFILE, "      o eqx line format:\n");
   fprintf(outFILE, "        - * = Match\n");
   fprintf(outFILE, "        - space = SNP/gap\n");
   fprintf(outFILE, "      o reference/query lines\n");
   fprintf(
      outFILE,
      "        - Without -print-positions\n"
   );
   fprintf(outFILE, "          - ID sequence\n");
   fprintf(
      outFILE,
      "        - With -print-positions\n"
   );

   fprintf(outFILE, "          - ID sequence end-base\n");

   if(defFormat == defFasta)
      fprintf(outFILE, "    -format-fasta: [Yes]\n");
   else fprintf(outFILE, "    -format-fasta: [No]\n");

   fprintf(
      outFILE, "      o Save alignment as a fasta file.\n"
   );
   fprintf(
      outFILE, "      o >id score first-base last-base\n"
   );
   fprintf(
      outFILE,"      o first-base is first aligned base\n"
   );
   fprintf(
      outFILE, "      o last-base is last aligned base\n"
   );

   /*****************************************************\
   * Fun-02 Sec-02 Sub-07:
   *  - Selecting alignment direction block
   \*****************************************************/

   #if defined SNPINSDEL
   #elif defined SNPDELINS
   #elif defined INSSNPDEL 
   #elif defined INSDELSNP
   #elif defined DELSNPINS
   #elif defined DELINSSNP
   #else

      fprintf(outFILE, "  Direction preferences:\n");
      if(defBestDir == defSnpInsDel)
         fprintf(outFILE, "    -snp-ins-del: [Yes]\n");
      else  fprintf(outFILE, "    -snp-ins-del: [No]\n");

      fprintf(
         outFILE,
         "      o For equal scores choose matches/SNPs"
      );
      fprintf(outFILE, " over\n        insertions and");
      fprintf(outFILE, " insertions over deletions.\n");

      if(defBestDir == defSnpDelIns)
         fprintf(outFILE, "    -snp-del-ins: [Yes]\n");
      else  fprintf(outFILE, "    -snp-del-ins: [No]\n");

      fprintf(
         outFILE,
         "      o For equal scores choose matches/SNPs"
      );
      fprintf(outFILE,"over\n        deletions and");
      fprintf(outFILE, " deletions over insertions.\n");

      if(defBestDir == defInsSnpDel)
         fprintf(outFILE, "    -ins-snp-del: [Yes]\n");
      else  fprintf(outFILE, "    -ins-snp-del: [No]\n");

      fprintf(
         outFILE,
         "      o For equal scores choose insertions over"
      );
      fprintf(outFILE, " matches/SNPs\n        and");
      fprintf(outFILE, " matches/SNPs over deletions.\n");

      if(defBestDir == defInsDelSnp)
         fprintf(outFILE, "    -ins-del-snp: [Yes]\n");
      else  fprintf(outFILE, "    -ins-del-snp: [No]\n");

      fprintf(
         outFILE,
         "      o For equal scores choose insertions over"
      );
      fprintf(outFILE, " deletions\n      and");
      fprintf(outFILE, " deletions over matches/SNPs.\n");

      if(defBestDir == defDelSnpIns)
         fprintf(outFILE, "    -del-snp-ins: [Yes]\n");
      else  fprintf(outFILE, "    -del-snp-ins: [No]\n");

      fprintf(
         outFILE,
         "      o For equal scores choose deletions over"
      );
      fprintf(outFILE, "matches/SNPs\n        and");
      fprintf(outFILE," matches/SNPs over insertions.\n");

      if(defBestDir == defDelInsSnp)
         fprintf(outFILE, "  -del-ins-snp: [Yes]\n");
      else  fprintf(outFILE, "  -del-ins-snp: [No]\n");

      fprintf(
         outFILE,
         "      o For equal scores choose deletions over"
      );
      fprintf(outFILE, " insertions\n        and");
     fprintf(outFILE, " insertions over matches/SNPs.\n");
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-03:
   ^  - Help message block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   helpOutputBlock:

     fprintf(outFILE, "  Help messages:\n");
     fprintf(outFILE, "    -h-all:\n");
     fprintf(
        outFILE,
        "      o Print the entire help message\n"
     );

   fprintf(
      outFILE,
      "    -v:\n      o Print alnSeq version.\n"
   );

   fprintf(
      outFILE,
      "    -flags:\n"
   );

   fprintf(
      outFILE,
      "      o Prints complier flags with descriptions\n"
   );

   fprintf(
      outFILE,
      "    -flags-only:\n"
   );

   fprintf(
      outFILE,
      "      o Prints the flags used and version number\n"
   );

   /*Output block*/
   fprintf(outFILE, "Output:\n");
   fprintf(
      outFILE,
      "  - Alignment to stdout or to file provided by"
   );
   fprintf(outFILE, " -out\n    and -alt-out.\n");

   return;
} /*printHelpMesg*/

/*-------------------------------------------------------\
| Fun-03: printCompilerSettings
|   - Prints out the compiler flags used
|   - Prints out what each possible compiler flag does
| Input:
|   - outFILE:
|     o File to print flags and flag descriptions to
|   - pDescBl:
|     o 1: print out descriptions for each flag (all)
|     o 0: Do not print out any flag descriptions
| Output:
|   - Prints flags and description to outFILE
\-------------------------------------------------------*/
void printCompileSettings(
   FILE *outFILE, /*Output file*/
   char pDescBl   /*not 0: print out flag descriptions*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC:
   '  - Print out compile flags used and what each
   '    compiler flag does.
   '  o fun-03 sec-01:
   '    - Print out the compiled settings
   '  o fun-03 sec-02:
   '    - Print out what each compiler flag does
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Print out the compiled settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(outFILE, "Compiler flags:\n");

   #if defined WORDS
      fprintf(outFILE, "   -DWORDS\n");
   #endif

   #if defined NOSEQCNVT
      fprintf(outFILE, "   -DNOSEQCNVT\n");
   #endif

   #if defined SNPINSDEL
      fprintf(outFILE, "   -DSNPINSDEL\n");
   #elif defined SNPDELINS
      fprintf(outFILE, "   -DSNPDELINS\n");
   #elif defined INSSNPDEL 
      fprintf(outFILE, "   -DINSSNPDEL\n");
   #elif defined INSDELSNP
      fprintf(outFILE, "   -DINSDELSNP\n");
   #elif defined DELSNPINS
      fprintf(outFILE, "   -DDELSNPINS\n");
   #elif defined DELINSSNP
      fprintf(outFILE, "   -DDELINSSNP\n");
   #endif

   if(!pDescBl) return;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Print out what each compiler flag does
   ^  o fun-03 sec-02 sub-01:
   ^    - Print out the header
   ^  o fun-03 sec-02 sub-02:
   ^    - Print out what the -DWORDS flag does
   ^  o fun-03 sec-02 sub-03:
   ^    - Print out what the -DBYTEMATRIX flag does
   ^  o fun-03 sec-02 sub-04:
   ^    - Print out what the -DHIRSCHTWOBIT flag does
   ^  o fun-03 sec-02 sub-05:
   ^    - Print out what the -DTWOBITMSW flag does
   ^  o fun-03 sec-02 sub-06:
   ^    - Print out what the -DNOGAPOPEN flag does
   ^  o fun-03 sec-02 sub-07:
   ^    - Print out what the -DSNPINSDEL flag does
   ^  o fun-03 sec-02 sub-08:
   ^    - Print out what the -DSNPDELINS flag does
   ^  o fun-03 sec-02 sub-09:
   ^    - Print out what the -DINSSNPDEL flag does
   ^  o fun-03 sec-02 sub-10:
   ^    - Print out what the -DINSDELSNP flag does
   ^  o fun-03 sec-02 sub-11:
   ^    - Print out what the -DDELSNPINS flag does
   ^  o fun-03 sec-02 sub-12:
   ^    - Print out what the -DDELINSSNP flag does
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - Print out the header
   \*****************************************************/

   fprintf(outFILE, "What each flag does:\n");

   fprintf(
      outFILE,
      "   - These flags are used to make alternative"
   );
   fprintf(outFILE, " versions\n     of alnSeq.\n");

   fprintf(
      outFILE,
      "   - Make alternate alnSeq:"
   );
   fprintf(
      outFILE,
      " make CFLAGS=\"flag1 flag2 ...\"\n"
   );

   /*****************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Print out what the -DWORDS flag does
   \*****************************************************/

   fprintf(outFILE, "   -DWORDS:\n");
   fprintf(
      outFILE,
      "     - Use a 128x128 (all ascii characters) matrix"
   );
   fprintf(
      outFILE,
      "\n       for the match and scoring matrix\n"
   );

   /*****************************************************\
   * Fun-03 Sec-02 Sub-03:
   *  - Print out what the -DNOSEQCNVT flag does
   \*****************************************************/

   fprintf(outFILE, "   -DNOSEQCNVT:\n");
   fprintf(
      outFILE,
      "     - Disables converting sequences to indexs\n"
   );
   fprintf(
      outFILE,
      "       This slows things down, but keeps case\n"
   );

   /*****************************************************\
   * Fun-03 Sec-02 Sub-07:
   *  - Print out what the -DSNPINSDEL flag does
   \*****************************************************/

   fprintf(outFILE, "   -DSNPINSDEL:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This"
   );
   fprintf(
      outFILE,
      " makes alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: SNPs/Matches, then insertions,"
   );
   fprintf(outFILE, "then\n       deletions.\n");

   /*****************************************************\
   * Fun-03 Sec-02 Sub-08:
   *  - Print out what the -DNSNPDELINS flag does
   \*****************************************************/

   fprintf(outFILE, "   -DSNPDELINS:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This"
   );
   fprintf(
      outFILE,
      " makes alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: SNPs/Matches, then deletinos, then"
   );
   fprintf(outFILE, "\n       insertions.\n");

   /*****************************************************\
   * Fun-03 Sec-02 Sub-09:
   *  - Print out what the -DINSSNPDEL flag does
   \*****************************************************/

   fprintf(outFILE, "   -DINSSNPDEL:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This"
   );
   fprintf(
      outFILE,
      " makes alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: insertions, then SNPs/Matches, then"
   );
   fprintf(outFILE, "\n       deletions.\n");

   /*****************************************************\
   * Fun-03 Sec-02 Sub-10:
   *  - Print out what the -DINSDELSNP flag does
   \*****************************************************/

   fprintf(outFILE, "   -DINSDELSNP:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This"
   );
   fprintf(
      outFILE,
      " makes alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: insertions, then deletions, then"
   );
   fprintf(outFILE, "\n       SNPs/Matches.\n");

   /*****************************************************\
   * Fun-03 Sec-02 Sub-11:
   *  - Print out what the -DDELSNPINS flag does
   \*****************************************************/

   fprintf(outFILE, "   -DDELSNPINS:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This"
   );
   fprintf(
      outFILE,
      " makes alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: deletions, then SNPs/Matches, then"
   );
   fprintf(outFILE, "\n       deletions.\n");

   /*****************************************************\
   * Fun-03 Sec-02 Sub-12:
   *  - Print out what the -DDELINSSNP flag does
   \*****************************************************/

   fprintf(outFILE, "   -DDELINSSNP:\n");
   fprintf(
      outFILE,
      "     - Disables the users ability to choose a "
   );
   fprintf(
      outFILE,
      "direction\n       when scores are equal. This"
   );
   fprintf(
      outFILE,
      " makes alignments\n       slightly faster.\n"
   );

   fprintf(
      outFILE,
      "     - Selects: deletions, then insertions, then"
   );
   fprintf(outFILE, "\n       SNPs/Matches.\n");

   return;
} /*printCompileSettings*/
