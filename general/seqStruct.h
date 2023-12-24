/*#########################################################
# Name: seqStruct
# Use:
#  o Holds functions for reading in or manipulating
#    sequences
# Libraries:
# C standard libraries:
#  - <stdlib.h>
#  - <stdio.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File:
'  o header:
'    - Has includes and default variables
'  o struct-01: seqStruct
'    - Holds sequence and length of a input sequence
'  o fun-01 addLineToBuffSeqFun:
'    - Add characters from file to buffer, if needed 
'      resize. This will only read in till the end of the
'      line
'  o fun-02 readFqSeq:
'    - Reads a fastq sequence from a fastq file
'  o fun-03 readFaSeq:
'    - Grabs the next read in the fasta file
'  o fun-04 complementBase:
'    - Returns the complement of a base
'  o fun-05 reverseComplementSeq:
'    - Reverse complement a sequence
'  o fun-06 reverseCStr;
'    - Reverse a c-string to be backwards (here for
'      Q-score entries)
'  o fun-07 initSeqST:
'     - Sets values in seqST to blank values
'  o fun-08 blankSeqST:
'    - Sets values in seqST to zero
'  o fun-09 freeSeqSTStack:
'    - Frees variables in an seqStruct (calls blankSeqST)
'  o fun-10 freeSeqST:
'    - Frees an seqST strucuter (calls fredSeqSTStack)
'  o fun-11 cpReadIdRPad:
'     - Copies read id to a buffer and adds in endIdC to
'       the end. If needed, this function will add right
'       padding of spaces to the end.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Header:
|  - Has includes and definede variables
\--------------------------------------------------------*/

#ifndef SEQSTRUCT_H
#define SEQSTRUCT_H

#include <stdio.h>
#include <stdlib.h>

/*--------------------------------------------------------\
| Struct-01: seqStruct
|  - Holds sequence and length of a input sequence
\--------------------------------------------------------*/
typedef struct seqStruct
{ /*seqStruct*/
  char *idCStr;               /*Id of th sequence*/
  unsigned long  lenIdUL;     /*Length of the sequence id*/
  unsigned long  lenIdBuffUL; /*Lenght of Id buffer*/

  char *seqCStr;              /*Sequence*/
  unsigned long  lenSeqUL;    /*Length of the sequence*/
  unsigned long  lenSeqBuffUL;/*Lenght of sequence buffer*/

  char *qCStr;                /*q-score entry*/
  unsigned long  lenQUL;      /*Length of the Q-score*/
  unsigned long  lenQBuffUL;  /*Length of Q-score buffer*/

  unsigned long  offsetUL;    /*Offset for an alignment*/
  unsigned long  endAlnUL;    /*Marks end of alignment*/
}seqStruct;

/*--------------------------------------------------------\
| Fun-01: addLineToBuffSeqFun
|  - Read line of characters into the buffer.If needed
|    this will resize the buffer.
| Input:
|  - buffCStr:
|    o Pionter to c-string to add buffer into.
|  - lenBuffUL:
|    o Pionter to the size of buffCStr
|  - curBuffUL:
|    o Pointer to the number of elements in buffCStr  
|  - resBuffUL:
|    o How much to expand buffer by when the buffer is
|      full
|  - inFILE:
|    o pointer to FILE handle to get the next lein from
| Output:
|   - Modifies:
|     o buffCStr to hold the next line.
|       - buffCStr is resizied if it is to small to hold
|         the next line.
|       - buffCStr + lenBuffUL - 2 will be '\0' or '\n'
|       - buffCStr set to 0 for memory allocation errors
|     o curBuffUL: Has the number of chars in the buffer
|     o lenBuffUL: Has the buffer size
|     o inFILE: Points to next line in file
|   - Returns:
|     o 0 if was end of file (EOF)
|     o 1 if read in the next line
|     o 64 if had a memory allocation error
\--------------------------------------------------------*/
static unsigned char addLineToBuffSeqFun(
    char **buffCStr,         /*Buffer to add data to*/
    unsigned long *lenBuffUL,/*Size of the buffer*/
    unsigned long *curBuffUL,/*Number of chars in buffer*/
    unsigned long resBuffUL, /*How much to resize buff by*/
    FILE *inFILE             /*File to grab data from*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: addLineToBuffSeqFun
   '  - Read line of characters into the buffer.If needed
   '    this will resize the buffer.
   '   o fun-01 sec-1:
   '     - variable declerations
   '   o fun-01 sec-2:
   '     - Check if need to resize the buffer
   '   o fun-01 sec-3:
   '     - Read in the next line in the buffer
   '   o fun-01 sec-4:
   '     - If at end of file, update read in lengths
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-1:
    ^  - variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    unsigned long spareBuffUL = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-2:
    ^  - Check if need to resize the buffer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*curBuffUL == 0 && *lenBuffUL > 0)
        spareBuffUL = *lenBuffUL;

    else if(*lenBuffUL == 0 || *curBuffUL -1 >= *lenBuffUL)
    { /*If need to resize the buffer (-1 for 0 index)*/
        *lenBuffUL += resBuffUL - 1;
           /*-1 to account for adding two index one items*/

        if(*lenBuffUL == 0)
          *buffCStr = malloc(sizeof(char) *(*lenBuffUL+1));
        else
        { /*Else I need to resize the buffer*/
            *buffCStr =
              realloc(
                *buffCStr,
                sizeof(char) * (*lenBuffUL + 1)
            ); /*Resize hte buffer*/
        } /*Else I need to resize the buffer*/
       
        if(*buffCStr == 0)
            return 64; /*Memory allocation error*/

        /*Amount of extra space in the buffer*/
        spareBuffUL = resBuffUL - 1; /*resBuff is index 1*/
    } /*If need to resize the buffer*/

    else
        spareBuffUL = *lenBuffUL - *curBuffUL;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-3:
    ^  - Read in the next line in the buffer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /*Set up marker to catch when read in the entire line*/
    *(*buffCStr + *lenBuffUL - 2) = '\0';
    tmpCStr = *buffCStr + *curBuffUL;

    while(fgets(tmpCStr, spareBuffUL, inFILE))
    { /*While I have lines to read*/

        if(*(*buffCStr + *lenBuffUL - 2) == '\0' ||
           *(*buffCStr + *lenBuffUL - 2) == '\n'
        ) { /*If read in the line*/

            /*Check if I used the entire buffer*/
            if(*(*buffCStr + *lenBuffUL - 2) == '\n')
                *curBuffUL = *lenBuffUL;
            else
            { /*Else only read in part of the buffer*/
                while(*tmpCStr != '\0')
                { /*While have characters in the buffer*/
                    ++tmpCStr;
                    ++(*curBuffUL);
                } /*While have characters in the buffer*/
            } /*Else only read in part of the buffer*/

            return 1; /*Finshed with the line*/
        } /*If read in the line*/

        /*Else resize the buffer*/
        *curBuffUL = *lenBuffUL - 1; /*Filled up buffer*/
            /*-1 to account for 0 index*/
        *lenBuffUL += resBuffUL - 1;
           /*-1 to account for adding two index 1 items*/

        *buffCStr =
          realloc(
            *buffCStr,
            sizeof(char) * (*lenBuffUL + 1)
        ); /*Resize the buffer*/

        if(*buffCStr == 0)
            return 64; /*Memory allocation error*/

        /*Reset my maker for entire line read in*/
        *(*buffCStr + *lenBuffUL - 2) = '\0';
        spareBuffUL = resBuffUL;
        tmpCStr = *buffCStr + *curBuffUL;
    } /*While I have lines to read*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-4:
    ^  - If at end of file, update read in lengths
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*(*buffCStr + *lenBuffUL - 2) == '\n')
        *curBuffUL = *lenBuffUL; /*used entire buffer*/
    else
    { /*Else only read in part of the buffer*/
        while(*tmpCStr != '\0')
        { /*While have characters in the buffer*/
            ++tmpCStr;
            ++(*curBuffUL);
        } /*While have characters in the buffer*/
    } /*Else only read in part of the buffer*/

    return 0; /*End of file*/
} /*addLineToBuffSeqFun*/

/*--------------------------------------------------------\
| Fun-02: readRefFqSeq
|  - Grabs the next read in the fastq file
| Input:
|  - fqFILE:
|    o Pointer to FILE handle to a fastq file to get the
|      next read from
|  - seqST:
|    o Pointer to seqStruct structure to hold the next read
| Output:
|  - Modifies:
|    o refStruct to hold the read in fastq entry & sets
|      its pointers
|  - Returns:
|     o 0: if EOF
|     o 1: if succeded
|     o 2: If file was not a fastq file
|     o 130: If file had an invalide entry
|       - This error uses to flags, first it uses 2 to
|         specify that it is not a fastq file
|       - 2nd it uses 128 to specifty that it is not an
|         blank file
|     o 64: If malloc failed to find memory
\--------------------------------------------------------*/
static unsigned char readFqSeq(
  FILE *fqFILE,           /*fastq file with sequence*/
  struct seqStruct *seqST /*Will hold one fastq entry*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: readRefFqSeq
   '  -  Grabs the next read in the fastq file
   '  o fun-02 sec-1:
   '    - Variable declarations
   '  o fun-02 sec-2:
   '    - Check if need to allocate memory for buffer
   '  o fun-02 sec-3:
   '    - Read in the first data
   '  o fun-02 sec-4:
   '    - If not at file end, see if have the full entry
   '  o fun-02 sec-5:
   '    -Read till end of file, check if valid fastq entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-02 Sec-1:
    ^  - Variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char tmpC = 'C';           /*initilize sequence loop*/

    unsigned short extraBuffUS = 1024;
    unsigned long tmpBuffUL = 0;

    /*Holds number of lines in sequence entry*/
    unsigned char numLinesUC = 0;

    unsigned char errUC = 0;
    char *oldIterCStr = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-02 Sec-2:
    ^  - Read in the header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(fqFILE == 0)
        return 2;    /*No file provided*/

    seqST->lenIdUL = 0;

    errUC =
        addLineToBuffSeqFun(
            &seqST->idCStr,     /*C-string to hold header*/
            &seqST->lenIdBuffUL,/*Length of header buff*/
            &seqST->lenIdUL,    /*Number bytes in buffer*/
            extraBuffUS,     /*Amount to increase buff by*/
            fqFILE        /*Fastq file to get header from*/
    ); /*Get the header (will resize as needed)*/

    /*have EOF (0) or memory allocation error (64)?*/
    if(!(errUC & 1)) return errUC;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-02 Sec-3:
    ^  - Read in the sequence & spacer
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    seqST->lenSeqUL = 0;
    /*need to set this up so the loop does not error out*/
    oldIterCStr = &tmpC;

    while(*oldIterCStr != '+')
    { /*While I have not reached the spacer entry*/
        /*So can get to this position later*/
        tmpBuffUL = seqST->lenSeqUL;

        errUC =
            addLineToBuffSeqFun(
                &seqST->seqCStr,    /*Buff to copy seq to*/
                &seqST->lenSeqBuffUL,/*Size of seq buffer*/
                &seqST->lenSeqUL,    /*Lenth of sequence*/
                extraBuffUS, /*How much to resize buff by*/
                fqFILE       /*Fastq file with sequence*/
        ); /*Get the header*/

        /*Set up new lines to be removed on the next read*/
        oldIterCStr = seqST->seqCStr + seqST->lenSeqUL - 1;

        while(*oldIterCStr < 33)
        { /*While removing new lines*/
            --seqST->lenSeqUL;
            --oldIterCStr;
        } /*While removing new lines*/

        /*Check for memory errors (64) & invalid entries*/
        if(errUC & 64) return errUC; 
        if(errUC == 0) return 2 + 128;

        /*Get on first character in the new buffer*/
        oldIterCStr = seqST->seqCStr + tmpBuffUL;
        ++numLinesUC; /*number of '\n' in sequence entry*/
    } /*While I have not reached the spacer entry*/

    --numLinesUC; /*Account for the overcounting*/

    oldIterCStr = seqST->seqCStr + seqST->lenSeqUL;

    while(*oldIterCStr != '+')
    { /*While not at start of spacer entry*/
        --oldIterCStr;
        --seqST->lenSeqUL; /*Acount for white space*/
    } /*While not at start of spacer entry*/

    *oldIterCStr = '\0';  /*Set spacer to null*/
    --oldIterCStr;
    --seqST->lenSeqUL;    /*Acount for white space*/

    while(*oldIterCStr < 33)
    { /*While have white space at end to remove*/
        *oldIterCStr = '\0';
        --oldIterCStr;
        --seqST->lenSeqUL; /*Acount for white space*/
    } /*While have white space at end to remove*/
    
    ++seqST->lenSeqUL; /*Account for being one base off*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-02 Sec-4:
    ^  - Read in the Q-score entry
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    seqST->lenQUL = 0;

    while(numLinesUC > 0)
    { /*While I need to read in the Q-score entry*/
        errUC =
            addLineToBuffSeqFun(
              &seqST->qCStr,   /*Buffer for q-score entry*/
              &seqST->lenQBuffUL,/*Size of buffer*/
              &seqST->lenQUL,   /*Length of Q-score entry*/
              extraBuffUS,      /*For reszing buffer*/
              fqFILE            /*Fastq file with q-score*/
        ); /*Get the header*/

        /*Set up white space to removed on next read*/
        oldIterCStr = seqST->qCStr + seqST->lenQUL - 1;

        while(*oldIterCStr < 33) /*removes white space*/
        { /*While removing new lines*/
            --seqST->lenQUL;
            --oldIterCStr;
        } /*While removing new lines*/

        /*Check for errors; memory=64/invalid file=130*/
        if(errUC & 64) return errUC;
        if(errUC == 0) return 2 + 128;

        --numLinesUC; /*Number of lines to read in*/
    } /*While I need to read in the Q-score entry*/

    /*Remove any white space at the end*/
    oldIterCStr = seqST->qCStr + seqST->lenQUL;

    while(*oldIterCStr < 33)
    { /*While have white space at end to remove*/
        *oldIterCStr = '\0';
        --oldIterCStr;
        --seqST->lenQUL; /*Acount for white space*/
    } /*While have white space at end to remove*/

    return errUC; /*Is 1 for another entry or 0 for EOF*/
} /*readFqSeq*/

/*--------------------------------------------------------\
| Fun-03 TOC: readFaSeq
|  -  Grabs the next read in the fasta file
| Input:
|  - faFILE:
|    o Pointer to FILE handle to a fasta file to get the
|      next read from
|  - seqST:
|    o Pointer to seqStruct structure to hold the next read
| Output:
|  - Modifies:
|    o seqST to hold one fasta entry
|  - Returns:
|     o 0: if EOF
|     o 1: if succeded
|     o 2: for an invalid fasta entry
|     o 64: If malloc failed to find memory
| Note:
|   - This will remove new lines from the sequence.
|   - This will only remove spaces or tabs at the end of
|     each sequence entry, so "atg atc \n" will got to
|     "atc atc".
\--------------------------------------------------------*/
static unsigned char readFaSeq(
  FILE *faFILE,           /*Fasta file with sequence*/
  struct seqStruct *seqST /*Will hold one fastq entry*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: readFaSeq
   '  -  Grabs the next read in the fasta file
   '  o fun-03 sec-1:
   '    - Variable declarations
   '  o fun-03 sec-2:
   '    - Check if need to allocate memory for buffer
   '  o fun-03 sec-3:
   '    - Read in the sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-03 Sec-1:
    ^  - Variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Used To initilize the sequence loop*/
    char tmpC = 'C';
    unsigned long tmpUL = 0;

    unsigned short extraBuffUS = 1024;
    unsigned char errUC = 0;
    char *oldIterCStr = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-03 Sec-2:
    ^  - Read in the header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(faFILE == 0)
        return 2;    /*No file provided*/

    seqST->lenIdUL = 0;

    errUC =
        addLineToBuffSeqFun(
            &seqST->idCStr,     /*Buffer to hold header*/
            &seqST->lenIdBuffUL,/*Length of buffer*/
            &seqST->lenIdUL,    /*length of header*/
            extraBuffUS,     /*How much to resize buff by*/
            faFILE             /*Fasta file with seq id*/
    ); /*Get the header (will resize as needed)*/

    /*Check for EOF (0) or memory allocation error (64)*/
    if(!(errUC & 1)) return errUC;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-03 Sec-3:
    ^  - Read in the sequence
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*need to set this up so the loop does not error out*/
    oldIterCStr = &tmpC;

    /*Overwrite the older sequence entry*/
    seqST->lenSeqUL = 0;

    if(seqST->seqCStr != 0) *seqST->seqCStr = '\0';

    while(*oldIterCStr != '>')
    { /*While I have not reached the spacer entry*/
        errUC =
            addLineToBuffSeqFun(
                &seqST->seqCStr,     /*Holds sequence*/
                &seqST->lenSeqBuffUL,/*size of seq buff*/
                &seqST->lenSeqUL,    /*Length of sequence*/
                extraBuffUS,
                faFILE         /*Fasta file with sequence*/
        ); /*Get the header*/

        /*Sets up to remove newlines on the next read in*/
        oldIterCStr = seqST->seqCStr + seqST->lenSeqUL - 1;

        while(*oldIterCStr < 33)/*Remove end white space*/
        { /*While removing new lines*/
            --seqST->lenSeqUL;
            --oldIterCStr;
        } /*While removing new lines*/

        /*Check for errors: memory (64), EOF (0)*/
        if(errUC & 64) return errUC;
        if(errUC == 0) break;

        /*Get on first character in the new buffer*/
        oldIterCStr = seqST->seqCStr + seqST->lenSeqUL;

        /*Peek ahead in the file*/
        tmpUL = fread(&tmpC, 1, 1, faFILE);
        if(tmpUL == 0) break;
        fseek(faFILE, -1, SEEK_CUR);
        if(tmpC == '>') break;
    } /*While I have not reached the spacer entry*/

    /*Check if is a valid fasta entry*/
    if(seqST->seqCStr == 0 || *seqST->seqCStr == '\0')
        return 2;

    /*Make sure the new line at the end is removed*/
    oldIterCStr = seqST->seqCStr + seqST->lenSeqUL;

    while(*oldIterCStr < 33)
    { /*While I have white space at the end*/
        *oldIterCStr = '\0';
        --oldIterCStr;
        --seqST->lenSeqUL;
    } /*While I have white space at the end*/

    ++seqST->lenSeqUL; /*Account for being one base off*/

    return errUC; /*Is 1 for another entry or 0 for EOF*/
} /*readFaSeq*/

/*--------------------------------------------------------\
| Fun-04: complementBase
|  - Returns the complement of a base
| Input:
|  - baseC:
|    o base to get the complement of
| Output:
|  - Returns
|    o the complement of the input base (0 if invalid base)
\--------------------------------------------------------*/
#define complementBase(baseC)({\
    char tmpC = 0;\
    switch((baseC) & (~32)) /*Make sure is upper case*/\
    { /*switch, reverse complement*/\
        case 'A': tmpC = 'T'; break;\
        case 'C': tmpC = 'G'; break;\
        case 'G': tmpC = 'C'; break;\
        case 'T': tmpC = 'A'; break;\
        case 'U': tmpC = 'A'; break;\
        case 'W': tmpC = 'W'; break;\
        case 'S': tmpC = 'S'; break;\
        case 'M': tmpC = 'K'; break;/*A/C to T/G*/\
        case 'K': tmpC = 'M'; break;/*T/G to A/C*/\
        case 'R': tmpC = 'Y'; break;/*A/G to T/C*/\
        case 'Y': tmpC = 'R'; break;/*T/C to A/G*/\
        case 'B': tmpC = 'V'; break;/*C/G/T to G/C/A*/\
        case 'D': tmpC = 'H'; break;/*G/T/A to C/A/T*/\
        case 'H': tmpC = 'D'; break;/*C/A/T to G/T/A*/\
        case 'V': tmpC = 'B'; break;/*A/C/G to T/G/C*/\
        case 'N': tmpC = 'N'; break;/*A/C/G/T to A/C/G/T*/\
        default: tmpC = 0;\
    } /*switch, reverse complement*/\
    \
    tmpC;\
}) /*complementBase*/

/*--------------------------------------------------------\
| Fun-05: reverseComplementSeq
|  - Reverse complement a sequence
| Input:
|  - seqST:
|    o Pointer to seqStruct structure with a sequence to
|      reverse complement
| Output:
|  - Modfies
|    o Reverse complements the sequence entry and if
|      reverses the Q-score entry
\--------------------------------------------------------*/
static void reverseComplementSeq(
  struct seqStruct *seqST  /*Sequence to reverse comp*/
){
    char *seqCStr = seqST->seqCStr;
    char *endCStr = seqST->seqCStr + seqST->lenSeqUL - 1;
    char *qCStr = 0;
    char *qEndCStr = 0;
    char swapC = 0;

    if(seqST->qCStr != 0 && *seqST->qCStr != '\0')
    { /*If have a Q-score entry*/
      qCStr = seqST->qCStr;
      qEndCStr = seqST->qCStr + seqST->lenQUL;
    } /*If have a Q-score entry*/

    while(endCStr > seqCStr)
    { /*While I have bases to reverse complement*/
        swapC = *seqCStr;
        *seqCStr = complementBase(*endCStr);
        *endCStr = complementBase(swapC);
        
        if(seqST->qCStr != 0 && *seqST->qCStr != '\0')
        { /*If I also need to swap Q-scores*/
          swapC = *qCStr;
          *qCStr = *qEndCStr;
          *qEndCStr = swapC;

          ++qCStr;
          --endCStr;
        } /*If I also need to swap Q-scores*/
        
        ++seqCStr;
        --endCStr;
    } /*While I have bases to reverse complement*/

    /*Check if ended on same base (if so complement base)*/
    if(endCStr == seqCStr)
        *seqCStr = complementBase(*seqCStr);

    return;
} /*reverseComplementSeq*/

/*--------------------------------------------------------\
| Fun-06: reverseCStr
|  - Reverse a c-string to be backwards
|    (here for Q-score entries)
| Input:
|  - inCStr:
|    o c-string to reverse
|  - lenCStrUI:
|    o lenght of inCStr (as index 1)
| Output:
|  - Modifies
|    o inCStr to be backwards (end at start, start at end)
\--------------------------------------------------------*/
static void reverseCStr(
  char *inCStr,       /*C-string to refeverse*/
  unsigned long lenCStrUI/*Length of input string index 1*/
){
    char *endCStr = inCStr + lenCStrUI - 1;
    char swapC = 0;

    while(endCStr > inCStr)
    { /*While I have bases to reverse complement*/
        swapC = *inCStr;

        *inCStr = *endCStr;
        *endCStr = swapC;
        ++inCStr;
        --endCStr;
    } /*While I have bases to reverse complement*/

    return;
} /*reverseCStr*/

/*--------------------------------------------------------\
| Fun-07: initSeqST
|  - Sets vlues in seqST to zero
| Input:
|  - seqSTPtr:
|    o Pointer to an seqStruct structuer to initialize
| Output:
|  - Modifies
|    o all values in seqST to be 0
\--------------------------------------------------------*/
#define initSeqST(seqSTPtr){\
   (seqSTPtr)->idCStr = 0;\
   (seqSTPtr)->lenIdUL = 0;\
   (seqSTPtr)->lenIdBuffUL = 0;\
   \
   (seqSTPtr)->seqCStr = 0;\
   (seqSTPtr)->lenSeqUL = 0;\
   (seqSTPtr)->lenSeqBuffUL = 0;\
   \
   (seqSTPtr)->qCStr = 0;\
   (seqSTPtr)->lenQUL = 0;\
   (seqSTPtr)->lenQBuffUL = 0;\
   \
   (seqSTPtr)->offsetUL = 0;\
   (seqSTPtr)->endAlnUL = 0;\
} /*initSeqST*/

/*--------------------------------------------------------\
| Fun-08: blankSeqST
|  - Sets values in seqST to blank values
| Input:
|  - seqSTPtr:
|    o Pointer to an seqStruct structuer to blank
| Output:
|  - Modifies
|    o Sets id, sequence, and Q-score entreis to start
|      with '\0' and the id, sequence, and Q-score lengths
|      to 0. This does not change the buffer lengths.
\--------------------------------------------------------*/
#define blankSeqST(seqSTPtr){\
   if((seqSTPtr)->idCStr != 0) *(seqSTPtr)->idCStr = '\0';\
   (seqSTPtr)->lenIdUL = 0;\
   \
   if((seqSTPtr)->seqCStr != 0) *(seqSTPtr)->seqCStr='\0';\
   (seqSTPtr)->lenSeqUL = 0;\
   \
   if((seqSTPtr)->qCStr != 0) *(seqSTPtr)->qCStr = '\0';\
   (seqSTPtr)->lenQUL = 0;\
   \
   (seqSTPtr)->offsetUL = 0;\
   (seqSTPtr)->endAlnUL = 0;\
} /*blankSeqST*/

/*--------------------------------------------------------\
| Fun-09: freeSeqSTStack
|  - Frees the variables in an seqStruct strucuter
| Input:
|  - seqSTPtr:
|    o Pointer to an seqStruct structuer that has variables
|      to free
| Output:
|  - Frees
|    o seqSTPtr-idCStr and sets to 0
|    o seqSTPtr-seqCStr and sets to 0
|    o seqSTPtr-qCStr and sets to 0
\--------------------------------------------------------*/
#define freeSeqSTStack(seqSTPtr){\
   if((seqSTPtr) != 0)\
   { /*If: I have something to free*/\
      blankSeqST((seqSTPtr));\
      \
      if((seqSTPtr)->idCStr != 0)\
         free((seqSTPtr)->idCStr);\
      if((seqSTPtr)->seqCStr != 0)\
         free((seqSTPtr)->seqCStr);\
      if((seqSTPtr)->qCStr != 0)\
         free((seqSTPtr)->qCStr);\
      \
      (seqSTPtr)->idCStr = 0;\
      (seqSTPtr)->seqCStr = 0;\
      (seqSTPtr)->qCStr = 0;\
   } /*If: I have something to free*/\
} /*freeSeqSTStack*/

/*--------------------------------------------------------\
| Fun-10: freeSeqST
|  - Frees the seqST strucuter
| Input:
|  - seqSTPtr:
|    o Pointer to an seqStruct structuer to free
| Output:
|  - Frees
|    o seqST from memory and sets to 0
\--------------------------------------------------------*/
#define freeSeqST(seqSTPtr){\
   freeSeqSTStack(seqSTPtr);\
   \
   if((seqSTPtr) != 0)\
   { /*If: I have something to free*/\
      free(seqSTPtr);\
      (seqSTPtr) = 0;\
   } /*If: I have something to free*/\
} /*freeSeqST*/

/*--------------------------------------------------------\
| Fun-11 TOC: cpReadIdRPad
|  - Copies read id to a buffer and adds in endIdC to
|    the end. If needed, this function will add right
|    padding of spaces to the end.
| Input:
|  - seqST:
|    o Pointer to seqStruct with read id to copy
|  - buffStr:
|    o c-string to copy read id to
|  - endIdC:
|    o character to add in to mark the end of the read id
|  - padRI:
|    o Amount of padding to add to the right of the read
|      id
| Output:
|  - Modifies:
|    o buffStr to hold the sequence id + endIdC. If the
|      id + endIdC is shorter than padRI, the copied id is
|      padded with spaces on the right till it is the same
|      size as padRI.
|  - Returns:
|    o Pointer to end of copied id or padding if padding
|      was applied.
\--------------------------------------------------------*/
static char * cpReadIdRPad(
   struct seqStruct *seqST, /*Has read id to copy*/
   char *buffStr,           /*Buffer to add read id to*/
   char endIdC,    /*Char to add to end of id (0 to skip)*/
   int padRI       /*Padding to add to right of id*/
){
   int idCntI = 0;
   char *tmpStr = seqST->idCStr;

   if(*tmpStr == '>') ++tmpStr; /*Get of header*/

   while(*tmpStr != '\0')
   { /*While I have a read id to copy over*/
      *buffStr = *tmpStr;
      ++buffStr;
      ++tmpStr;
      ++idCntI;
   } /*While I have a read id to copy over*/

   if(endIdC)
   { /*If adding in a character after the read id*/
      *buffStr = endIdC;
      ++buffStr; /*If not a null*/
      ++idCntI;
   } /*If adding in a character after the read id*/

   while(idCntI < padRI)
   { /*While I have padding to add to the end*/
      *buffStr = ' ';
      ++buffStr;
      ++idCntI;
   } /*While I have padding to add to the end*/

   return buffStr;
} /*cpReadIdRPad*/

#endif
