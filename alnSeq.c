/*#########################################################
# Name: alignSeq
# Use:
#  - Runs a Needlman Wunsch or Smith Waterman alignment on
#    a pair of fasta files
# Includes:
#  - "waterman.h"
#  - "needleman.h"
#  - "sequenceFun.h"
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
#  - <string.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <stdint.h>
#########################################################*/

#include "waterman.h"
#include "needleman.h"
#include "sequenceFun.h"
#include <string.h>

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
    struct alnSet *alnSetST // Aligment settings
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '    fun-01 sec-1: Variable declerations
   '    fun-01 sec-2: Look through user input
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
   '    - Make aligned query sequence
   '  o main sec-08:
   '    - Make aligned reference sequence
   '  o main sec-09:
   '    - Print out the single alignment and clean up
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

   // Holds the reference sequence
   struct seqStruct refST;
   struct seqStruct queryST;

   // Temporay variables for getting the sequence from the
   // fasta file (holds sequence buffer size)
   uint32_t lenSeqBuffUI = 0;

   // Output aligned sequences
   char *queryAlnCStr = 0;
   char *refAlnCStr = 0;
   long bestScoreL = 0;

   // Caputures error type from functions
   unsigned char errUC = 0;
    
   // For holding settings
   struct alnSet alnSetST;

   // For holding alignment output
   struct alnMatrixStruct *alnMtrxST = 0;
   struct alnStruct *alnST = 0;

   FILE *faFILE = 0;
   FILE *outFILE = 0;

   char *helpCStr = "\
       \n alnSeq -query query.fasta -ref ref.fasta options\
       \n Use:\
       \n  - Runs a Needleman-Wunsch or Smith Waterman\
       \n    alignment on input sequences\
       \n Input:\
       \n  -query: [Required]\
       \n    o Fasta sequence to align.\
       \n  -ref: [Required]\
       \n    o Second fasta sequence to align.\
       \n  -out: [stdout]\
       \n    o File to output alignment\
       \n  -gapopen: [-1]\
       \n    o Cost of starting an indel (as integer).\
       \n    o A negative value is a penalty, while\
       \n      postives values are correct.\
       \n    o Posivitve values favor gaps, while negative\
       \n      values favors SNPs.\
       \n  -gapextend: [-4]\
       \n    o Cost of extending an indel.\
       \n    o Similar to gapopen.\
       \n  -score-matrix: [ENDNAFULL]\
       \n    o File with matrix to use. It should have one\
       \n      line for each possible score. Scores not\
       \n      supplied will be set to defaults.\
       \n      (see scoring-matrix.txt file in source).\
       \n    o Format:\
       \n      a t -4\
       \n      a a 5 \
       \n      \\\\ This is a comment\
       \n  -use-water: [Needleman Wunsch]\
       \n     o Use a Waterman Smith alignment instead of\
       \n       the Needleman Wunsch\
       \n     o This only tracks one best alignment, not\
       \n       all best (or all) alignments.\
       \n  -line-wrap: [59]\
       \n    o Maximum characters per line in the output\
       \n      alignment file.\
       \n    o The minimum line wrap is 42.\
       \n  -multi-base-water: [No]\
       \n    o For a Waterman Smith alignment, print out\
       \n      the best alignment for each base.\
       \n  -min-score: [100]\
       \n   - Minimum score needed to keep an non-best\
       \n     alignment when -multi-base-water is used.\
       \n  -min-bases: [50]\
       \n   - Minimum number unique bases needed for\
       \n     reference or query to keep an non-best\
       \n     alignment when -multi-base-water is used.\
       \n  -enable-match-priority: [No]\
       \n     o Always select matches, even when an indel\
       \n       would give a better score.\
       \n     o This does not apply to SNPs.\
       \n  -match-ins-del: [match-ins-del]\
       \n     o For equal scores choose matches/SNPs over\
       \n       insertions and insertions over deletions.\
       \n  -match-del-ins: [match-del-ins]\
       \n     o For equal scores choose matches/SNPs over\
       \n       deletions and deletions over insertions.\
       \n  -ins-match-del: [ins-match-del]\
       \n     o For equal scores choose insertions over\
       \n       matches/SNPs and matchs/SNPs over dels.\
       \n  -del-match-ins: [del-match-ins]\
       \n     o For equal scores choose deletions over\
       \n       matches/SNPs and matchs/SNPs over ins.\
       \n  -ins-del-match: [ins-del-match]\
       \n     o For equal scores choose insertions over\
       \n       deletions and deletions over matchs/SNPs.\
       \n  -del-ins-match: [del-ins-match]\
       \n     o For equal scores choose deletions over\
       \n       insertions and ins over matches/SNPs.\
       \n Output:\
       \n  - Alignment with query on top, reference, and\
       \n    error line to stdout or specified file.\
       \n  - Error line format:\
       \n    o I = Insertion\
       \n    o D = Deletion\
       \n    o X = mismatch\
       \n    o = = match\
       \n    o S is soft maksing for query and reference\
       \n    o s is soft maksing for query only\
       \n    o P is soft maksing for reference only\
       ";

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-2:
   ^  - Read in user input and check input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&alnSetST);

   inputCStr =
       checkInput(
           &lenArgsInt,
           argsCStr,
           &refFileCStr,
           &queryFileCStr,
           &outFileCStr,
           &alnSetST
   ); // Get the user input

   if(inputCStr != 0)
   { // If had problematic input
        if(strcmp(inputCStr, "-h") == 0 ||
           strcmp(inputCStr, "--h") == 0 ||
           strcmp(inputCStr, "-help") == 0 ||
           strcmp(inputCStr, "--help") == 0 ||
           strcmp(inputCStr, "help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpCStr);
            exit(0);
        } /*If user wanted the help message*/

        if(strcmp(inputCStr, "-V") == 0 ||
           strcmp(inputCStr, "-v") == 0 ||
           strcmp(inputCStr, "--V") == 0 ||
           strcmp(inputCStr, "--v") == 0 ||
           strcmp(inputCStr, "--version") == 0 ||
           strcmp(inputCStr, "--Version") == 0 ||
           strcmp(inputCStr, "-version") == 0 ||
           strcmp(inputCStr, "-Version") == 0
        ) { /*if the user wanted the version number*/
            fprintf(
                stdout,
                "alignSeq from findCoInft version: %.8f\n",
                defVersion
            ); /*Print out the closest thing to a version*/
            exit(0);
        } /*Else if the user wanted the version number*/

        // If file with scoring matrix was invalid
        else if(strcmp(inputCStr, "-score-matrix") == 0)
            exit(1);

        else if(inputCStr != 0)
        { /*If user had invalid input*/
            fprintf(
                stderr,
                "%s\n%s is invalid\n",
                helpCStr,
                inputCStr
            ); /*Print out the problem*/
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
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   faFILE = fopen(refFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       printf(
         "Reference (-ref %s) could not be opend\n",
         refFileCStr
       );

       exit(-1);
   } // If reference file could not be opened

   lenSeqBuffUI = 0;

   errUC = 
       readFaSeq(
           faFILE,
           &refST.idCStr,
           &refST.lenIdUC,
           &refST.seqCStr,
           &lenSeqBuffUI,
           &refST.lenSeqUI
   ); // Read in the reference sequence

   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       freeSeqST(
         refST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // refST on stack (do not free)
       );

       printf(
         "Reference (-ref %s) is not valid\n",
         refFileCStr
       );

       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
       freeSeqST(
         refST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // refST on stack (do not free)
       );

       printf("Memory allocation error\n");
       exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-04:
   ^  - read in the query sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   faFILE = fopen(queryFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       freeSeqST(
         refST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // refST on stack (do not free)
       );

       printf(
         "Query (-query %s) could not be opend\n",
         queryFileCStr
       );

       exit(-1);
   } // If reference file could not be opened

   lenSeqBuffUI = 0;

   errUC = 
       readFaSeq(
           faFILE,
           &queryST.idCStr,
           &querySt.lenIdUC,
           &queryST.seqCStr,
           &lenSeqBuffUI,
           &queryST.lenSeqUI
   ); // Read in the reference sequence

   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       freeSeqST(
         refST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // refST on stack (do not free)
       );

       freeSeqST(
         queryST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // queryST on stack (do not free)
       );

       printf(
         "Query (-query %s) is not valid\n",
         refFileCStr
       );
       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
      freeSeqST(
         refST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // refST on stack (do not free)
       );

       freeSeqST(
         queryST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // queryST on stack (do not free)
       );

      printf("Memory allocation error\n");
      exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-05:
   ^  - Do the alingment
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   if(alnSetST.useNeedleBl != 0)
     alnMtrxST = NeedlemanAln(queryST, refST, &alnSetST);

   else alnMtrxST = WatermanAln(queryST,refST,&alnSetST);

   if(alnMtrxST == 0)
   { // If did not have enough memory
      freeSeqST(
         refST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // refST on stack (do not free)
       );

       freeSeqST(
         queryST,
         1, // Free the sequence array
         1, // Free the read id array
         0  // queryST on stack (do not free)
       );

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
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   if(outFileCStr != 0) outFILE = fopen(outFileCStr, "w");
   else outFILE = stdout;

   if(
     alnSetST->mulitBaseWaterBl == 1 &&
     alnSetST.useNeedleBl == 0
   ) { // If doing a multi-aligment extraction
     if(outFileCStr == 0)
       outFileCStr = "out";

     printAlnWaterAlns(
       alnMtrxST,
       &queryST,
       &refST,
       alnSetST,
       outFileCStr // Prefix to name everything
     );

     freeAlnMatrixST(alnMtrxST); // No longer need
     goto alnSeqCleanup;
   } // If doing a multi-aligment extraction
   
   alnST =
     cnvtDirMatrixToAlnAry(
       &refST,
       &queryST,
       alnMtrxST->dirMatrixST,
       alnMtrxST->bestScoreST,
       1                          // Include soft masking
   );

   bestScoreL = alnSt->bestScoreST->scoreL;
   freeAlnMatrixST(alnMtrxST); // No longer need

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-07:
   ^  - Make aligned query sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   queryAlnCStr = cnvtAlnToSeq(queryST, 1, alnST);

   if(queryAlnCStr == 0)
   { // If had a memory allocation error
     alnAryToLetter(refST->seqCStr,queryST->seqCStr,alnST);
     fprintf(outFILE, "%s\n", alnErrUCAry);

     fprintf(
       stderr,
       "Memory error, while outputing sequences\n"
     );

     fprintf(stderr, "Printed alignment array\n");

     freeSeqST(refST, 1, 1, 0);
     freeSeqST(queryST, 1, 1, 0);
     freeAlnST(alnST);

     fclose(outFILE);
     outFILE = 0;

     exit(-1);
   } // If had a memory allocation error

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-08:
   ^  - Make aligned reference sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   refAlnCStr = cnvtAlnToSeq(refST, 0, alnST);

   if(refAlnCStr == 0)
   { // If had a memory allocation error
     alnAryToLetter(refST->seqCStr,queryST->seqCStr,alnST);
     fprintf(outFILE, "%s\n%s\n", queryAlnCStr, alnErrUCAry);

     fprintf(
          stderr,
          "Memory error, while outputing sequences\n"
     );

     fprintf(stderr, "Printed alignment array\n");

     freeSeqST(refST, 1, 1, 0);
     freeSeqST(queryST, 1, 1, 0);
     freeAlnST(alnST);

     fclose(outFILE);
     outFILE = 0;

     exit(-1);
   } // If had a memory allocation error

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-08:
   ^  - Print out the alignment and clean up
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   alnAryToLetter(refST->seqCStr, queryST->seqCStr, alnST);

   printAln(
     outFILE,
     queryST->idCStr,
     queryAlnCStr,
     refST->idCStr,
     refAlnCStr,
     scoreL,
     alnSetST->lineWrapUS,
     alnST
  );
     
   alnSeqCleanup:

   fclose(outFILE);
   outFILE = 0;

   free(queryAlnCStr);
   free(refAlnCStr);
   freeAlnST(alnST); // NEED TO SET UP
   freeSeqST(refST, 1, 1, 0);
   freeSeqST(queryST, 1, 1, 0);

   exit(0);
} // main

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
    struct alnSet *alnSetST // Aligment settings
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

    for(int intArg = 1; intArg < *lenArgsInt; intArg++)
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
            alnSetST->gapStartPenaltyI =
                strtol(singleArgCStr, &tmpCStr, 10);

        else if(strcmp(tmpCStr, "-gapextend") == 0)
            alnSetST->gapExtendPenaltyI =
                strtol(singleArgCStr, &tmpCStr, 10);

        else if(strcmp(tmpCStr, "-line-wrap") == 0)
            cStrToUSht(singleArgCStr,alnSetST->lineWrapUS);

        else if(strcmp(tmpCStr, "-use-water") == 0)
        { // Else if disabling match priority
            alnSetST->useNeedleBl = !defUseNeedle;
            --intArg;
        } // Else if disabling match priority

        else if(strcmp(tmpCStr, "-multi-base-water") == 0)
        { // Else if doing more than the best alignment
          alnSetST->mulitBaseWaterBl = !defMultiBaseWater;
          alnSetST->useNeedleBl = !defUseNeedle;
          --intArg;
        } // Else if doing more than the best alignment

        else if(strcmp(tmpCStr, "-min-score") == 0)
          cStrToUInt(singleArgCStr, alnSetST->minScoreUI);

        else if(strcmp(tmpCStr, "-min-bases") == 0)
          cStrToUInt(singleArgCStr, alnSetST->minBasesUI);

        else if(
          strcmp(tmpCStr, "-enable-match-priority") == 0
        ) { // Else if disabling match priority
            alnSetST->matchPriorityBl = !defMatchPriority;
            --intArg;
        } // Else if disabling match priority

        else if(strcmp(tmpCStr, "-match-ins-del") == 0)
        { // Else if wants matches->insertions->deletions
            alnSetST->diagnolPriorityC = 0;
            alnSetST->topPriorityC = 1;
            alnSetST->leftPriorityC = 2;
            --intArg;
        } // Else if wants matches->insertions->deletions

        else if(strcmp(tmpCStr, "-match-del-ins") == 0)
        { // Else if wants matches->deletions->insertions
            alnSetST->diagnolPriorityC = 0;
            alnSetST->topPriorityC = 2;
            alnSetST->leftPriorityC = 1;
            --intArg;
        } // Else if wants matches->deletions->insertions

        else if(strcmp(tmpCStr, "-ins-match-del") == 0)
        { // Else if wants insertions->matches->deletions
            alnSetST->diagnolPriorityC = 1;
            alnSetST->topPriorityC = 0;
            alnSetST->leftPriorityC = 2;
            --intArg;
        } // Else if wants insertions->matches->deletions

        else if(strcmp(tmpCStr, "-del-match-ins") == 0)
        { // Else if wants deletions->matches->insertions
            alnSetST->diagnolPriorityC = 1;
            alnSetST->topPriorityC = 2;
            alnSetST->leftPriorityC = 0;
            --intArg;
        } // Else if wants deletions->matches->insertions

        else if(strcmp(tmpCStr, "-ins-del-match") == 0)
        { // Else if wants insertions->deletions->matches
            alnSetST->diagnolPriorityC = 2;
            alnSetST->topPriorityC = 0;
            alnSetST->leftPriorityC = 1;
            --intArg;
        } // Else if wants insertions->deletions->matches

        else if(strcmp(tmpCStr, "-del-ins-match") == 0)
        { // Else if wants deletions->insertions->matches
            alnSetST->diagnolPriorityC = 2;
            alnSetST->topPriorityC = 1;
            alnSetST->leftPriorityC = 0;
            --intArg;
        } // Else if wants deletions->insertions->matches

        else if(strcmp(tmpCStr, "-score-matrix") == 0)
        { // else if the user supplied a scoring matrix
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
              readInScoreFile(alnSetST, scoreFILE);

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

        // Move to the parameter, so next input is a flag
        intArg++;
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/
