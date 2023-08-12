/*#########################################################
# Name: alignSeq
# Use:
#  - Runs a Needlman Wunsch or Smith Waterman alignment on
#    a pair of fasta files
# Includes:
#  - "hirschberg.h"
#  - "waterman.h"
#  o "needleman.h"
#  o "generalAlnFun.h"
#  o "alnStruct.h"
#  o "alnMatrixStruct.h"
#  o "twoBitArrays.h"
#  o "scoresST.h"
#  o "seqStruct.h"
#  o "alnSetStruct.h"
#  o "alnSeqDefaults.h"
#  o "twoBitArrays.h"
# C standard libraries:
#  o <string.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <stdint.h>
#########################################################*/

#include "hirschberg.h"
#include "waterman.h"  // hirschberg.h
//#include "needleman.h" // hirschberg.h
//#include <string.h> // in waterman.h

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOP: Start Of Program
'  - main: Run this program
'  - fun-01 checkInput;
'    o Gets user input
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output: Modifies: Each input variable to hold user input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to "-score-matrix" for invalid scoring
|       matrix input
|  - Prints to stdout when the scoring file is invalid
\--------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,        // Number arguments user input
    char *argsCStr[],       // Array with user arguments
    char **refFileCStr,     // file name of reference file
    char **queryFileCStr,   // File name of the query file
    char **outFileCStr,     // Name of the output file
    char **scoreMtrxFileStr, /*Holds scoring matrix file*/
    struct alnSet *settings // Aligment settings
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '    fun-01 sec-1: Variable declerations
   '    fun-01 sec-2: Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Prints:
|    o help message to outFILE.
\--------------------------------------------------------*/
void printHelpMesg(
   FILE *outFILE,
   char breifBl   /*Print a shorter help message*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printHelpMesg
   '  - Prints out the help message to outFILE
   '  o fun-03 Sec-01:
   '    - Usage block
   '  o fun-03 Sec-02:
   '    - Input block
   '  o fun-03 Sec-03:
   '    - Output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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

   // User input
   char *refFileCStr = 0;
   char *queryFileCStr = 0;
   char *outFileCStr = 0;
   char *inputCStr = 0;
   char *scoreMtrxFileStr = 0;

   // Holds the reference sequence
   struct seqStruct refST;
   struct seqStruct queryST;

   // Output aligned sequences
   char *queryAlnCStr = 0;
   char *refAlnCStr = 0;

   // Caputures error type from functions
   unsigned char errUC = 0;
   long bestScoreL = 0;  // Allows me to free score matrix

   // For holding settings
   struct alnSet settings;

   // For holding alignment output
   struct alnMatrixStruct *alnMtrxST = 0;
   struct alnStruct *alnST = 0;

   FILE *faFILE = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-2:
   ^  - Read in user input and check input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&settings);

   inputCStr =
       checkInput(
           &lenArgsInt,
           argsCStr,
           &refFileCStr,
           &queryFileCStr,
           &outFileCStr,
           &scoreMtrxFileStr,
           &settings
   ); // Get the user input

   if(inputCStr != 0)
   { // If had problematic input
        if(strcmp(inputCStr, "-h") == 0 ||
           strcmp(inputCStr, "--h") == 0 ||
           strcmp(inputCStr, "-help") == 0 ||
           strcmp(inputCStr, "--help") == 0 ||
           strcmp(inputCStr, "help") == 0
        ) { /*If user wanted the help message*/
            printHelpMesg(stdout, 1);/*Breif help message*/
            exit(0);
        } /*If user wanted the help message*/

        if(strcmp(inputCStr, "-h-all") == 0 ||
           strcmp(inputCStr, "--h-all") == 0 ||
           strcmp(inputCStr, "-help-all") == 0 ||
           strcmp(inputCStr, "--help-all") == 0 ||
           strcmp(inputCStr, "help-all") == 0
        ) { /*If user wanted the help message*/
            printHelpMesg(stdout, 0); /*full help message*/
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
            ); /*Print out the closest thing to a version*/
            exit(0);
        } /*Else if the user wanted the version number*/

        // If file with scoring matrix was invalid
        else if(strcmp(inputCStr, "-score-matrix") == 0)
            exit(1);

        else if(inputCStr != 0)
        { /*If user had invalid input*/
            printHelpMesg(stderr, 1); /*short help*/
            fprintf(stderr, "%s is invalid\n", inputCStr);
            exit(1); /*Let user know their was an error*/
        } /*If user had invalid input*/
   } // If had problematic input

   if(outFileCStr != 0)
   { // If printing output to a file
        outFILE = fopen(outFileCStr, "w");

        if(outFILE == 0)
        { // If an invalid output file
            printf(
              "Output (-out %s) file is invalid.\n",
              outFileCStr
            ); // Let user know about the invalid file

            exit(-1);
        } // If an invalid output file

        fclose(outFILE);
        outFILE = 0;
   } // If printing output to a file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-3:
   ^  - read in the reference sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   faFILE = fopen(refFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       fprintf(
         stderr,
         "Reference (-ref %s) could not be opend\n",
         refFileCStr
       );

       exit(-1);
   } // If reference file could not be opened

   // Read in the reference sequence
   initSeqST(&refST);
   errUC = readFaSeq(faFILE, &refST);
   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       freeSeqST(&refST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Reference (-ref %s) is not valid\n",
         refFileCStr
       );

       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
       freeSeqST(&refST, 0); // 0 to makr on the stack
       fprintf(stderr, "Memory allocation error\n");
       exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-04:
   ^  - read in the query sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   faFILE = fopen(queryFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       freeSeqST(&refST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Query (-query %s) could not be opend\n",
         queryFileCStr
       );

       exit(-1);
   } // If reference file could not be opened


   // Read in the query sequence
   initSeqST(&queryST);
   errUC = readFaSeq(faFILE, &queryST);
   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       freeSeqST(&refST, 0); // 0 to makr on the stack
       freeSeqST(&queryST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Query (-query %s) is not valid\n",
         refFileCStr
       );
       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
      freeSeqST(&refST, 0); // 0 to makr on the stack
      freeSeqST(&queryST, 0); // 0 to makr on the stack
      fprintf(stderr, "Memory allocation error\n");
      exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-05:
   ^  - Do the alingment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Right know these are hardoced in, but at some piont
   // it might be nice to allow the user the manipulate.
   // I would need to set up the Waterman and Needleman
   // alignments to handle this

   seqToLookupIndex(&refST);
   seqToLookupIndex(&queryST);

   queryST.endAlnUL = queryST.lenSeqUL - 1;
   queryST.offsetUL = 0;

   refST.endAlnUL = refST.lenSeqUL - 1;
   refST.offsetUL = 0;

   if(outFileCStr != 0) outFILE = fopen(outFileCStr, "w");

   else
   { // Else putting output to stdout
     outFILE = stdout;
     outFileCStr = "out";
   } // Else putting output to stdout

   if(settings.useNeedleBl != 0)
     alnMtrxST = NeedlemanAln(&queryST, &refST, &settings);

   else if(settings.useWaterBl != 0)
     alnMtrxST =
       WatermanAln(&queryST,&refST,&settings,outFileCStr);

   else if(settings.useHirschBl != 0)
   { /*Else if doing an Hirschberg alignment*/
     alnST = Hirschberg(&refST, &queryST, &settings);
     if(alnST == 0) goto alignmentFailed;

     /* The Hirschberg returns an alignment structure,
     `  instead of an directional matrix.
     */
     if(settings.useHirschBl != 0) goto noDirMatrix;
   } /*Else if doing an Hirschberg alignment*/

   if(alnMtrxST == 0)
   { // If did not have enough memory
      alignmentFailed:

      freeSeqST(&refST, 0); // 0 to makr on the stack
      freeSeqST(&queryST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Memory allocation error in aligment step\n"
       );
       exit(-1);
   } // If did not have enough memory

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-06:
   ^  - Print out multi-alignments or find the single
   ^    alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
     settings.multiBaseWaterBl == 1 &&
     settings.matrixScanBl == 0
   ) printAltWaterAlns(
        alnMtrxST,
        &queryST,
        &refST,
        &settings,
        outFileCStr /*Prefix to name files*/
      );
   
   alnST =
     dirMatrixToAlnST(
       &refST,
       &queryST,
       alnMtrxST->dirMatrixST,
       &alnMtrxST->bestScoreST
   );

   bestScoreL = alnMtrxST->bestScoreST.scoreL;
   freeAlnMatrixST(alnMtrxST); // No longer need
   alnMtrxST = 0;

   if(alnST == 0)
   { // IF i falied to make an alignment array
      freeSeqST(&refST, 0); // 0 to makr on the stack
      freeSeqST(&queryST, 0); // 0 to makr on the stack

       fprintf(
         stderr,
         "Memory error when making alignment array\n"
       );

       exit(-1);
   } // IF i falied to make an alignment array

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-07:
   ^  - Print out the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   noDirMatrix: /*For aligners that return an alnStruct*/

   lookupIndexToSeq(&refST);
   lookupIndexToSeq(&queryST);

   if(
      printAln(
         outFILE,
         outFileCStr,
         &refST,
         &queryST,
         alnST,
         bestScoreL,
         &settings,
         scoreMtrxFileStr /*File name for scoring matrix*/
      )
   ){ /*If could not print out the alignmnet*/
      fprintf(stderr, "Failed to print alignment\n");

      fclose(outFILE);
      free(queryAlnCStr);
      free(refAlnCStr);
      freeAlnST(alnST, 1); /*NEED TO SET UP*/

      freeSeqST(&refST, 0);   /*0 to specify on stack*/
      freeSeqST(&queryST, 0); /* 0 to do a stack free*/

      exit(1);
   } /*If could not print out the alignmnet*/

   fclose(outFILE);
   free(queryAlnCStr);
   free(refAlnCStr);
   freeAlnST(alnST, 1);

   freeSeqST(&refST, 0);   /*0 to specify on stack*/
   freeSeqST(&queryST, 0); /* 0 to do a stack free*/

   exit(0);
} /*main*/

/*--------------------------------------------------------\
| Output: Modifies: Each input variable to hold user input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to "-score-matrix" for invalid scoring
|       matrix input
|  - Prints to stdout when the scoring file is invalid
\--------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,        // Number arguments user input
    char *argsCStr[],       // Array with user arguments
    char **refFileCStr,     // file name of reference file
    char **queryFileCStr,   // File name of the query file
    char **outFileCStr,     // Name of the output file
    char **scoreMtrxFileStr, /*Holds scoring matrix file*/
    struct alnSet *settings // Aligment settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '    fun-01 sec-1: Variable declerations
   '    fun-01 sec-2: Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    char *singleArgCStr = 0;
    unsigned long scoreFileErrUL = 0;
    FILE *scoreFILE = 0;  // For loading the scoring matrix

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-01 Sec-2: Look through user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg += 2)
    { /*loop through all user input arguments*/
        // The 0 index holds the program name
        singleArgCStr = *(argsCStr +intArg + 1);// argument
        tmpCStr = *(argsCStr + intArg);         // Paramter

        if(strcmp(tmpCStr, "-ref") == 0)
            *refFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-query") == 0)
            *queryFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-out") == 0)
            *outFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-gapopen") == 0)
            settings->gapOpenI =
                strtol(singleArgCStr, &tmpCStr, 10);

        else if(strcmp(tmpCStr, "-gapextend") == 0)
          settings->gapExtendI =
                strtol(singleArgCStr, &tmpCStr, 10);

        else if(strcmp(tmpCStr, "-print-aligned") == 0)
        { /*Else if printing only the aligned region*/
           settings->pFullAlnBl = 0;
           --intArg;
        } /*Else if printing only the aligned region*/

        else if(strcmp(tmpCStr, "-print-unaligned") == 0)
        { /*Else if printing only the aligned region*/
           settings->pFullAlnBl = 1;
           --intArg;
        } /*Else if printing only the aligned region*/

        else if(strcmp(tmpCStr, "-print-positions") == 0)
        { /*Else if printing out the base positions*/
           settings->pBasePosBl = 1;
           --intArg;
        } /*Else if printing out the base positions*/

        else if(strcmp(tmpCStr, "-no-positions") == 0)
        { /*Else if not printing out the base positions*/
           settings->pBasePosBl = 0;
           --intArg;
        } /*Else if not printing out the base positions*/

        else if(strcmp(tmpCStr, "-line-wrap") == 0)
          cStrToUSht(singleArgCStr, &settings->lineWrapUS);

        else if(strcmp(tmpCStr, "-use-needle") == 0)
        { // Else if disabling match priority
            settings->useNeedleBl = 1;
            settings->useWaterBl = 0;
            settings->useHirschBl = 0;
            --intArg;
        } // Else if disabling match priority

        else if(strcmp(tmpCStr, "-use-water") == 0)
        { // Else if doing a waterman smith alignment
            settings->useNeedleBl = 0;
            settings->useWaterBl = 1;
            settings->useHirschBl = 0;
            --intArg;
        } // Else if doing a waterman smith alignment

        else if(strcmp(tmpCStr, "-use-hirschberg") == 0)
        { // Else if doing a Hirshberg alignment
            settings->useNeedleBl = 0;
            settings->useWaterBl = 0;
            settings->useHirschBl = 1;
            --intArg;
        } // Else if doing a Hirshberg alignment

        else if(strcmp(tmpCStr, "-format-expand-cig") == 0)
        { /*Else if using the expanded cigar format*/
           settings->formatFlag = defExpandCig;
           --intArg;
        } /*Else if using an expanded cigar format*/

        else if(strcmp(tmpCStr, "-format-emboss") == 0)
        { /*Else if using an EMBOSS like format*/
           settings->formatFlag = defEMBOSS;
           --intArg;
        } /*Else if using an EMBOSS like format*/

        else if(strcmp(tmpCStr, "-format-clustal") == 0)
        { /*Else if using clustal format*/
           settings->formatFlag = defClustal;
           --intArg;
        } /*Else if using clustal format*/

        else if(strcmp(tmpCStr, "-format-fasta") == 0)
        { /*Else if outputing fasta format*/
           settings->formatFlag = defFasta;
           --intArg;
        } /*Else if outputing fasta format*/

        else if(
          strcmp(tmpCStr, "-query-ref-scan-water") == 0
        )
        { // Else if doing more than the best alignment
          settings->multiBaseWaterBl = 1;
          settings->refQueryScanBl = 1;
          settings->matrixScanBl = 0;

          settings->useNeedleBl = 0;
          settings->useWaterBl = 1;
          settings->useHirschBl = 0;

          --intArg;
        } // Else if doing more than the best alignment

        else if(strcmp(tmpCStr, "-matrix-scan-water") == 0)
        { // Else if doing a matrix scan
          settings->multiBaseWaterBl = 1;
          settings->refQueryScanBl = 0;
          settings->matrixScanBl = 1;

          settings->useNeedleBl = 0;
          settings->useWaterBl = 1;
          settings->useHirschBl = 0;

          --intArg;
        } // Else if doing a matrix scan

        else if(strcmp(tmpCStr, "-min-score") == 0)
          cStrToUInt(singleArgCStr, &settings->minScoreUI);

        // Not used
        //else if(strcmp(tmpCStr, "-min-bases") == 0)
        //  cStrToUInt(singleArgCStr, &settings->minBasesUI);

        else if(strcmp(tmpCStr, "-match-ins-del") == 0)
        { // Else if wants matches->insertions->deletions
            settings->diagnolPriorityC = 0;
            settings->topPriorityC = 1;
            settings->leftPriorityC = 2;
            --intArg;
        } // Else if wants matches->insertions->deletions

        else if(strcmp(tmpCStr, "-match-del-ins") == 0)
        { // Else if wants matches->deletions->insertions
            settings->diagnolPriorityC = 0;
            settings->topPriorityC = 2;
            settings->leftPriorityC = 1;
            --intArg;
        } // Else if wants matches->deletions->insertions

        else if(strcmp(tmpCStr, "-ins-match-del") == 0)
        { // Else if wants insertions->matches->deletions
            settings->diagnolPriorityC = 1;
            settings->topPriorityC = 0;
            settings->leftPriorityC = 2;
            --intArg;
        } // Else if wants insertions->matches->deletions

        else if(strcmp(tmpCStr, "-del-match-ins") == 0)
        { // Else if wants deletions->matches->insertions
            settings->diagnolPriorityC = 1;
            settings->topPriorityC = 2;
            settings->leftPriorityC = 0;
            --intArg;
        } // Else if wants deletions->matches->insertions

        else if(strcmp(tmpCStr, "-ins-del-match") == 0)
        { // Else if wants insertions->deletions->matches
            settings->diagnolPriorityC = 2;
            settings->topPriorityC = 0;
            settings->leftPriorityC = 1;
            --intArg;
        } // Else if wants insertions->deletions->matches

        else if(strcmp(tmpCStr, "-del-ins-match") == 0)
        { // Else if wants deletions->insertions->matches
            settings->diagnolPriorityC = 2;
            settings->topPriorityC = 1;
            settings->leftPriorityC = 0;
            --intArg;
        } // Else if wants deletions->insertions->matches

        else if(strcmp(tmpCStr, "-score-matrix") == 0)
        { // else if the user supplied a scoring matrix
            *scoreMtrxFileStr = singleArgCStr;
            scoreFILE = fopen(singleArgCStr, "r");

            if(scoreFILE == 0)
            { // If I could not open the scoring file
                return tmpCStr; // So user knows invalid

                fprintf(stderr,
                  "-score-matrix %s is not an file\n",
                  singleArgCStr
                ); /*Print out the problem*/
            } // If I could not open the scoring file

            scoreFileErrUL =
              readInScoreFile(settings, scoreFILE);

            if(scoreFileErrUL != 0)
            { // If the scoring file had an error
              fprintf(
                stderr,
                "Invalid line (%lu) in -score-matrix %s\n",
                scoreFileErrUL,
                singleArgCStr
              ); /*Print out the problem*/

                return tmpCStr;    // invalid file
            } // If the scoring file had an error
        } // else if the user supplied a scoring matrix
            
        else return tmpCStr; // Invalid parameter
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/

/*--------------------------------------------------------\
| Output:
|  - Prints:
|    o help message to outFILE.
\--------------------------------------------------------*/
void printHelpMesg(
   FILE *outFILE,
   char breifBl   /*Print a shorter help message*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: printHelpMesg
   '  - Prints out the help message to outFILE
   '  o fun-03 Sec-01:
   '    - Usage block
   '  o fun-03 Sec-02:
   '    - Input block
   '  o fun-03 Sec-03:
   '    - Output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Usage block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      outFILE,
      "alnSeq -query query.fasta -ref ref.fasta"
   );
   fprintf(outFILE, " [options...]\n");

   fprintf(outFILE, "Use:\n");
   fprintf(
     outFILE,
     "  - Does a pairwise alignment on two input sequences"
   );
   fprintf(outFILE, "\n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Input block
   ^  o fun-03 sec-02 sub-01:
   ^     - General input/output block (includes line wrap)
   ^  o fun-03 sec-02 sub-02:
   ^     - Alignment algorithim selection block
   ^  o fun-03 sec-02 sub-03:
   ^     - Alignment paramaters block
   ^  o fun-03 sec-02 sub-04:
   ^    - File output settings (non-format)
   ^  o fun-03 sec-02 sub-05:
   ^     - File output format block
   ^  o fun-03 sec-02 sub-06:
   ^     - Waterman specific paramters block
   ^  o fun-03 sec-02 sub-07:
   ^     - Selecting alignment direction block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - General input/output block (includes line wrap)
   \******************************************************/

   fprintf(outFILE, "Input:\n");

   fprintf(outFILE, "  -query: [Required]\n");
   fprintf(
      outFILE,
      "    o Fasta file with query sequence.\n"
   );

   fprintf(outFILE, "  -ref: [Required]\n");
   fprintf(
      outFILE,
      "    o Fasta file with reference sequence.\n"
   );

   fprintf(outFILE, "  -out: [stdout]\n");
   fprintf(outFILE,"    o File to output alignment to.\n");

   /******************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Alignment algorithim selection block
   \******************************************************/

   if(defUseNeedle)
      fprintf(outFILE, "  -use-needle: [Yes]\n");
   else fprintf(outFILE, "  -use-needle: [No]\n");

   fprintf(
      outFILE,
      "    o Do a Needleman Wunsch alignment.\n" 
   );

   if(defUseWater)
      fprintf(outFILE, "  -use-water: [Yes]\n");
   else fprintf(outFILE, "  -use-water: [No]\n");

   fprintf(
      outFILE,
      "    o Do a Waterman Smith alignment.\n"
   );

   if(defUseHirsch)
       fprintf(outFILE, "  -use-hirschberg: [Yes]\n");
   else fprintf(outFILE, "  -use-hirschberg: [No]\n");

   fprintf(outFILE, "    o Do a Hirschberg alignment.\n");

   /******************************************************\
   * Fun-03 Sec-02 Sub-03:
   *  - Alignment paramaters block
   \******************************************************/

   fprintf(outFILE, "  -gapopen: [%i]\n", defGapOpen);
   fprintf(
      outFILE,
      "    o Cost of starting an indel (as integer). A "
   );
  fprintf(outFILE,"negative\n      value is a penalty.\n");

   fprintf(outFILE, "  -gapextend: [%i]\n", defGapExtend);
   fprintf(
      outFILE,
      "    o Cost of extending an indel one base (< 0 is"
   );
   fprintf(outFILE, " penalty).\n");

   fprintf(
      outFILE,
      "  -score-matrix: [%s]\n",
      defMatrixNameStr
   );
   fprintf(
      outFILE,
      "    o File with scoring matrix to use. It should"
   );
   fprintf(
      outFILE,
      " have one\n      line for each score to change.\n"
   );
   fprintf(
      outFILE,
      "    o (see scoring-matrix.txt for an example).\n"
   );

   if(breifBl) goto helpOutputBlock;

   /******************************************************\
   * Fun-03 Sec-02 Sub-04:
   *  - File output settings (non-format)
   \******************************************************/

   if(defPAln)
      fprintf(outFILE,"  -print-aligned: [Yes]\n");
   else fprintf(outFILE,"  -print-aligned: [No]\n");

   fprintf(
      outFILE,
      "    o Only print the aligned regions of each"
   );
   fprintf(outFILE, " sequence.\n");

   if(!defPAln)
      fprintf(outFILE,"  -print-unaligned: [Yes]\n");
   else fprintf(outFILE,"  -print-unaligned: [No]\n");

   fprintf(
      outFILE,
      "    o Print out the entire reference and query"
   );
   fprintf(outFILE, " sequence.\n");

   if(defPPos)
      fprintf(outFILE, "  -print-positions: [Yes]\n");
   else fprintf(outFILE, "  -print-positions: [No]\n");

   fprintf(
      outFILE,
      "    o Print the starting and ending position for"
   );
   fprintf(outFILE, " each line\n");
   fprintf(outFILE, "      in the alignment.\n");
   fprintf(
      outFILE,
     "    o Clustal format only prints the ending position"
   );
   fprintf(outFILE, ".\n");

   if(!defPPos)
      fprintf(outFILE, "  -no-positions: [Yes]\n");
   else fprintf(outFILE, "  -no-positions: [No]\n");

   fprintf(
      outFILE,
      "    o Turns off -print-positions\n"
   );
   fprintf(
      outFILE,
      "    o EMBOSS format always prints positions.\n"
   );


   fprintf(outFILE, "  -line-wrap: [%i]\n", defLineWrap);
   fprintf(
      outFILE,
      "    o Maximum characters per line in output file.\n"
   );
   fprintf(outFILE, "    o Input 0 for no line wrap.\n");
   fprintf(
      outFILE,
      "    o Minimum line wrap is 10 for fasta (header not"
   );
   fprintf(outFILE, "\n");
   fprintf(
      outFILE,
      "      wrapped), 32 for clustal, & 42 for expanded"
   );
   fprintf(outFILE, " cigar\n      and EMBOSS.\n");

   /******************************************************\
   * Fun-03 Sec-02 Sub-05:
   *  - File output format block
   \******************************************************/

   if(defFormat == defExpandCig)
      fprintf(outFILE, "  -format-expand-cig: [Yes]\n");
   else fprintf(outFILE, "  -format-expand-cig: [No]\n");

   fprintf(
      outFILE,
      "    o Prints the reference sequence, then query"
   );
   fprintf(outFILE, " sequence,\n");
   fprintf(outFILE, "      and then the eqx line.\n");
   fprintf(outFILE, "    o eqx line format:\n");
   fprintf(outFILE, "      - I = Insertion\n");
   fprintf(outFILE, "      - D = Deletion\n");
   fprintf(outFILE, "      - X = mismatch\n");
   fprintf(outFILE, "      - = = match\n");
   fprintf(outFILE, "      - S = soft mask\n");
   fprintf(outFILE, "    o reference/query lines\n");
   fprintf(
      outFILE,
      "      - Without -print-positions\n"
   );

   fprintf(outFILE, "        - Ref: sequence\n");
   fprintf(outFILE, "        - Qry: sequence\n");
   fprintf(
      outFILE,
      "      - With -print-positions\n"
   );
   fprintf(
      outFILE,
      "        - Ref: start-base sequence end-base\n"
   );
   fprintf(
      outFILE,
      "        - Qry: start-base sequence end-base\n"
   );


   if(defFormat == defEMBOSS)
      fprintf(outFILE, "  -format-emboss: [Yes]\n");
   else fprintf(outFILE, "  -format-emboss: [No]\n");

   fprintf(
      outFILE,
      "    o Prints the reference sequence, then eqx line,"
   );
   fprintf(outFILE,"and then\n      the query sequence.");
   fprintf(outFILE, "\n    o eqx line format:\n");
   fprintf(outFILE, "      - | = Match\n");
   fprintf(outFILE, "      - space = SNP/gap\n");
   fprintf(outFILE, "    o reference/query lines\n");
   fprintf(
      outFILE,
      "      - ID: start-base sequence end-base\n"
   );

   if(defFormat == defClustal)
      fprintf(outFILE, "  -format-clustal: [Yes]\n");
   else fprintf(outFILE, "  -format-clustal: [No]\n");

   fprintf(
      outFILE,
      "    o Prints the reference sequence, then query"
   );
   fprintf(outFILE, " sequence,\n");
   fprintf(outFILE, "      and then the eqx line.\n");
   fprintf(outFILE, "    o eqx line format:\n");
   fprintf(outFILE, "      - * = Match\n");
   fprintf(outFILE, "      - space = SNP/gap\n");
   fprintf(outFILE, "    o reference/query lines\n");
   fprintf(
      outFILE,
      "      - Without -print-positions\n"
   );
   fprintf(outFILE, "        - ID sequence\n");
   fprintf(
      outFILE,
      "      - With -print-positions\n"
   );

   fprintf(outFILE, "        - ID sequence end-base\n");

   if(defFormat == defFasta)
      fprintf(outFILE, "  -format-fasta: [Yes]\n");
   else fprintf(outFILE, "  -format-fasta: [No]\n");

   fprintf(
      outFILE, "    o Save alignment as a fasta file.\n"
   );
   fprintf(
      outFILE, "    o >id score first-base last-base\n"
   );
   fprintf(
      outFILE, "    o first-base is first aligned base\n"
   );
   fprintf(
      outFILE, "    o last-base is last aligned base\n"
   );

   /******************************************************\
   * Fun-03 Sec-02 Sub-06:
   *  - Waterman specific paramters block
   \******************************************************/

   if(defQueryRefScan)
       fprintf(outFILE,"  -query-ref-scan-water: [Yes]\n");
   else fprintf(outFILE,"  -query-ref-scan-water: [No]\n");

   fprintf(outFILE, "    o Waterman alignment only.\n");
   fprintf(
      outFILE,
      "    o prints out the best score for each"
   );
   fprintf(outFILE, " reference and\n      query base.\n");

   if(defMatrixScan)
      fprintf(outFILE, "  -matrix-scan-water: [Yes]\n");
   else fprintf(outFILE, "  -matrix-scan-water: [No]\n");

   fprintf(outFILE, "    o Waterman alignment only.\n");
   fprintf(
      outFILE,
      "    o Print out every full aligment above a min"
   );
   fprintf(
      outFILE,
      " score\n      in a dirction matrix.\n");

   fprintf(outFILE, "  -min-score: [%i]\n", defMinScore);
   fprintf(outFILE, "    o Waterman alignment only.\n");
   fprintf(
      outFILE,
      "    o Minimum score needed to keep an alternative\n"
   );
   fprintf(outFILE, "      alignment.\n");

   /******************************************************\
   * Fun-03 Sec-02 Sub-07:
   *  - Selecting alignment direction block
   \******************************************************/

   if(
      defTopPriority == 1 &&
      defDiagnolPriority == 0 &&
      defLeftPriority == 2
   ) fprintf(outFILE, "  -match-ins-del: [Yes]\n");

   else  fprintf(outFILE, "  -match-ins-del: [No]\n");

   fprintf(
      outFILE,
      "    o For equal scores choose matches/SNPs over "
   );
   fprintf(outFILE, "insertions\n      and");
   fprintf(outFILE, " insertions over deletions.\n");

   if(
      defTopPriority == 2 &&
      defDiagnolPriority == 0 &&
      defLeftPriority == 1
   ) fprintf(outFILE, "  -match-del-ins: [Yes]\n");

   else  fprintf(outFILE, "  -match-del-ins: [No]\n");

   fprintf(
      outFILE,
      "    o For equal scores choose matches/SNPs over "
   );
   fprintf(outFILE, "deletions\n      and");
   fprintf(outFILE, " deletions over insertions.\n");

   if(
      defTopPriority == 0 &&
      defDiagnolPriority == 1 &&
      defLeftPriority == 2
   ) fprintf(outFILE, "  -ins-match-del: [Yes]\n");

   else  fprintf(outFILE, "  -ins-match-del: [No]\n");

   fprintf(
      outFILE,
      "    o For equal scores choose insertions over "
   );
   fprintf(outFILE, "matches/SNPs\n      and");
   fprintf(outFILE, " matches/SNPs over deletions.\n");

   if(
      defTopPriority == 2 &&
      defDiagnolPriority == 1 &&
      defLeftPriority == 0
   ) fprintf(outFILE, "  -del-match-ins: [Yes]\n");

   else  fprintf(outFILE, "  -del-match-ins: [No]\n");

   fprintf(
      outFILE,
      "    o For equal scores choose deletions over "
   );
   fprintf(outFILE, "matches/SNPs\n      and");
   fprintf(outFILE, " matches/SNPs over insertions.\n");

   if(
      defTopPriority == 0 &&
      defDiagnolPriority == 2 &&
      defLeftPriority == 1
   ) fprintf(outFILE, "  -ins-del-match: [Yes]\n");

   else  fprintf(outFILE, "  -ins-del-match: [No]\n");

   fprintf(
      outFILE,
      "    o For equal scores choose insertions over "
   );
   fprintf(outFILE, "deletions\n      and");
   fprintf(outFILE, " deletions over matches/SNPs.\n");

   if(
      defTopPriority == 1 &&
      defDiagnolPriority == 2 &&
      defLeftPriority == 0
   ) fprintf(outFILE, "  -del-ins-match: [Yes]\n");

   else  fprintf(outFILE, "  -del-ins-match: [No]\n");

   fprintf(
      outFILE,
      "    o For equal scores choose deletions over "
   );
   fprintf(outFILE, "insertions\n      and");
   fprintf(outFILE, " insertions over matches/SNPs.\n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Output block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   helpOutputBlock:

   if(breifBl)
     fprintf(
       outFILE,
       "  -h-all:\n      o Print the entire help message\n"
      );

   /*Output block*/
   fprintf(outFILE, "Output:\n");
   fprintf(
      outFILE,
      "  - Alignment to stdout or to file provided by"
   );
   fprintf(outFILE, " -out.\n");

   return;
} /*printHelpMesg*/
