# Use:

This is here to give you an idea of the functions and data
  structures you will need to use to use this code in your
  own code.

This is not the best guide and does not cover everything,
  but I hope it will help.

## Reading in sequences

Sequences are stored in the seqStruct, which is struct-01
  of seqStruct.h. seqStruct stores the sequence id,
  sequence, and q-score entries. It also stores the length
  of the id, sequence, and q-score and the size of each
  buffer storing the id, sequence, and q-score.

An seqStruct should be initialized with the initSeqST
  function (fun-08 in seqStruct.h) only for the first use.
  After initialization use the blankSeqST() function 
  (fun-10 seqStruct.h) to reset the structure to no
  sequence.

Use the freeSeqST function (fun-07 seqStruct.h) to clean
  up. freeSeqST takes a pointer to an seqStruct structure
  and a variable to mark if the structure is on the heap
  or stack. If heapBl is set to 1, then the input structure
  will be freed, but if not, then only the buffers will be
  freed. freeSeqST(seqStruct, 1) does not set the pointer
  to the seqStruct to 0, you must do this.

A sequence can be read in from a fastq file using the
  readFqSeq function (fun-03 in seqStruct.h) or fasta file
  using the readFaSeq function (fun-04 in seqStruct.h).
  Both functions take a pointer to the FILE with the
  sequence and a pointer to an seqStruct structure. The
  FILE pointer will be set to the next entry in the
  fasta/fastq file.

You can reverse complement a sequence using the
  reverseComplementSeq() function (fun-01 in seqStruct.h).

```
Here is an example.

FILE faFILE = fopen("sequence.fasta", "r");
struct seqStruct *sequence = 0;

if(faFILE == 0) return 0;

sequence=malloc(sizeof(struct seqStruct));

if(sequence == 0)
{
  fclose(faFILE);
  return 0;
}

initSeqST(&sequence);

switch(readFaSeq(faFILE, &sequence))
{ // Switch: Check what kind of error
  case 0:
    fclose(faFILE);
    freeSeqST(sequence, 1);
    sequence = 0;
    return 0;        // EOF

  case 1:
    fclose(faFILE);
    return sequence;

  case 2:
    fclose(faFILE);
    freeSeqST(sequence, 1);
    sequence = 0;
    return 0;        // invalid file

  case 64:
    fclose(faFILE);
    freeSeqST(sequence, 1);
    sequence = 0;
    return 0;       // memory error
} // Switch: Check what kind of error
```

## Alignment settings

Settings for the alignment settings are stored in the
  alnSet structure (struct-01 alnSetStruct.h). This
  structure holds the output line wrap size, which
  alignment to run, if you are doing a matrix scan, if
  you are doing a reference/query scan, the scoring matrix,
  the gap extension penalty, the gap starting penalty, the
  min score to keep an alternative alignment, and the order
  to prefer a single direction. The default values for
  these settings can all be modified in alnSeqDefaults.h.
  
You can initialize an alnSet structure using the
  initAlnSet() function (fun-01 in alnSetStruct.h). This 
  will set all values in the alnSet structure to the
  default values in alnSeqDefaults.h.

For changing the scoring matrix you can use the
  setBasePairScore() function (fun-02 in alnSetStruct.h) to
  change the value for individual base pairs in a matrix or
  use the readInScoreFile() function
  (fun-04 alnSetStruct.h) to read in a scoring matrix. This
  scoring matrix is different than the normal format, so
  see scoring-matrix.txt for an example of my format.

Changing the non-scoring values for an alnSet structure
  requires you to manually change them. If the value ends
  in Bl, then it should only ever be a 1 or 0.

For the matrix scan (matrixScanBl) and reference/query
  scan, make sure that multiBaseWaterBl is set to 1 and
  useNeedleBl is set to 0.

Finally you can free the alnSet structure using the 
  freeAlnSet() function (fun-03 in alnSetStruct.h). This
  takes in a pointer to an alnSet structure and a 1 or 0.
  If you use a 1, then the alnSet structure will not be
  freed, but if you use a 0 it will be freed. Just remember
  to set the pointer to null calling freeAlnSet(alnSet, 0);

```
An example

struct alnSet *settings = 0;
FILE *matrixFILE = fopen("scoring-matrix.txt", "r");

if(matrixFILE == 0) return ERROR;

settings = malloc(sizeof(struct alnSet));

if(settings == 0)
{
  fclose(matrixFILE);
  return ERROR;
}

initAlnSet(settings);
readInScoreFile(alnSetST, matrixFILE);

// Do something

fclose(matrixFILE);
freeAlnSet(settings, 0);

return 0;
```
  
## Doing an alignment

alnSeq can do either an Needleman Wunsch or an Waterman
  Smith alignment. The Needleman Wunsch alignment is done
  by the NeedlemanAln function (fun-01 in needle.h), while
  the Waterman Smith alignment is done by the WatermanAln
  function (fun-01 water.h). Both functions take an pointer
  to an seqStruct structure with a query sequence, an
  pointer to an seqStruct structure with a reference
  sequence, an pointer to a alignment settings structure.
  In addition the WatermanAln also takes in a prefix to
  name the output matrix scan file if a matrix scan is
  requested.

Both the Needleman Wunsch and Waterman Smith alignments
  will return an alnMatrixStruct structure
  (struct-01 alnMatrixStruct.h), which will contain the
  directional matrix as a two bit array, best score/ending
  score, and the extra scores kept for the
  reference/query scan (if was called for).

Once you are finished with the matrix you can free it with
  the freeAlnMatrixST() function (fun-02 in
  alnMatrixStruct.h). Just remember to set it to null after
  freeing.

```
struct seqStruct querySeq;
struct seqStruct refSeq;
struct alnSet settings;

struct alnMatrixStruct *matrix = 0;

FILE *faFILE = 0;

initAlnSet(&settings);
initSeqST(&querySeq);
initSeqST(&refSeq);

faFILE = fopen("query.fasta", "r");
if(faFILE == 0) return 0;
fclose(faFILE);

switch(readFaSeq(faFILE, &querySeq))
{
// Check to see if read in sequence & handle errors here
}

faFILE = fopen("ref.fasta", "r");

if(faFILE == 0)
{
// Handle file error here
}

fclose(faFILE);

switch(readFaSeq(faFILE, &refSeq))
{
// Check to see if read in sequence & handle errors here
}

// This is to avoid errors. Originally I wanted to be able
// to align in the middle of sequences, but I never set
// that option up fully. It would crash the current code.
queryST.endAlnUI = queryST.lenSeqUI;
queryST.offsetUI = 0;

refST.endAlnUI = refST.lenSeqUI;
refST.offsetUI = 0;

matrix = WatermanAln(querySeq, refSeqSt, settings, "");
// or matrix = NeedlemanAln(querySeq, refSeqSt, settings);

// Do something with matrix here

freeAlnMatrixST(alignmentMatrix);
freeAlnSet(settings, 0);
```
  
## Printing alignments

For printing alignments you have two options, the first
  is to print out the score, ending query base, ending
  reference base, backwards cigar, starting query base,
  and starting reference base. The second is to print out
  an alignment.

### Printing a cigar

You can print out the cigar entry using the printMatrixCig
  function (fun-03 water.h). This function takes in a 
  FILE to output to, a directional matrix in the form of
  a two bit array, the length of the reference, and the
  score of the alignment. The matrix should have a pointer
  to the index of the starting base, which can be set by
  using
  `moveXElmFromStart(alnMatrixStruct->dirMatrixST, index)`
  (Fun-15 in twoBitArrays.c). For the best alignment use
  alnMatrixStruct->bestScoreSt->indexUL.

The downside of printMatrixCig() is that it is not very
  well coded, so it will run a bit slow.

```
Printing out the best alignment as a cigar would look like

// Set up and run alignment and set up the output file

// Get the length of the aligned region of the reference
referenceLength = refSeq->enalnUI - refSeq->offset;

moveXElmFromStart(
  matrix->dirMatrixST,
  matrix->bestScoreST->indexUL
); // Move to the best scoring base pair in the matrix

printMatrixCig(
  outputFILE,
  matrix->dirMatrixST,
  referenceLength,
  matrix->bestScoreST->scoreL
); // Print out the backwards cigar
```

### Printing an alignment

To print out the alignment you first have to create the
  alignment array, the aligned reference sequence, the
  aligned query sequence, and convert the alignment array
  into human readable format. After that you can print out the
  alignment using the printAln function (fun-04 in
  alnStruct.h).

Converting the directional matrix to an alignment array 
  is done with the cnvtDirMatrixToAlnAry() function
  (fun-03 in alnStruct.h). This function takes in a pointer
  to an seqStruct structure with the reference sequence,
  an pointer to an seqStruct structure with the query
  sequence, an two bit direction matrix
  (alnMatrixStruct->dirMatrixST), a scoreStruct structure
  (alnMatrixStruct->bestScoreST), and a 1 (add soft masking
  to the ends) or 0.

After building the alignment array you can then free the
  alnMatrixStruct structure. Just remember to save the best
  score before freeing the alnMatrixStruct.

The aligned query and reference sequences can be built
  using the cnvtAlnAryToSeq() function (fun-01
  alnStruct.h). This takes in seqStruct structure with the
  sequence to make into an aligned sequence, a 1 (query) or
  0 (reference) to tell if working on the query or
  reference sequence, and an pointer to the alignment array
  returned form cnvtDirMatrixToAlnAry.

After aligning the reference and query sequences you can
  convert the alignment array to human readable format with
  the alnAryToLetter() function (fun-02 alnStruct.h). This
  takes a pointer to an seqStruct structure with the 
  reference sequence, an pointer to an seqStruct to the
  query sequence, and the alignment array structure from
  cnvtAlnAryToSeq.

At this point you are read to print out the alignment usin
  the printAln() function. Just make sure to free the
  alignment array structure (alnStruct) with
  freeAlnSt(alnStruct, 1);. Were 1 tells the function to
  free the structure and variables in the structure. A 0
  would just free the variable inside the structure.

```
char *alignedQuery = 0;
char *alignedRef = 0;
long bestScoreL = 0;
struct alnStruct *alignment = 0;

// Other variables, such as the output file

// Get the alignment array
alignment =
  cnvtDirMatrixToAlnAry(
    refSeq,
    querySeq,
    matrix->dirMatrixST,
    matrix->besScoreST,
    1
); // Build the alignment array 

bestScoreL->matrix->bestScoreST->scoreL.

freeAlnMatrixST(matrix);
matrix=0;

// Get the aligned query sequence
alignedQuery = cnvtAlnAryToSeq(querySeq, 1, alignment);

if(alignedQuery == 0)
{
  // Handle errors here
}

// Get the aligned reference sequence
alignedRef = cnvtAlnAryToSeq(querySeq, 1, alignment);

if(alignedRef == 0)
{
  // Handle errors here
}

// Printing out the structure
alnAryToLetter(refSeq, querySeq, alignment);

printAln(
  outFILE,
  querySeq->idCStr,     // c-string from seqStruct
  alignedQuery,         // c-string
  refSeq->idCStr,       // c-string form seqStruct
  alignedRef,           // c-string
  bestScoreL,           // Long
  settings->lineWrapUS, // unsigned short from alnSetStruct
  alignment             // alnStruct
  );

fclose(outFILE);

freeSeqST(querySeq);
freeSeqST(refSeq);

free(alignedQuery);
free(alignedRef);

freeAlnST(alignment, 1);
freeAlnSet(settings, 0);

return 0;
```

## Two bit arrays

The two bit array is an array of uint8_t's that each hold
  four two bit elements. This array is stored in a
  twoBitAry structure (struct-01 twoBitArrays.h). The
  twoBitAry structure has a pointer to the first unit8_t
  (limb) in the array (firstLimbUCPtr), a pointer to the 
  limb currently one (limbOnUCPtr), a counter telling which
  element we are on in the limb, and the length of the
  array.

You can create a twoBitAry structure using the
  makeTwoBitArry function (fun-12 in twoBitArrys.h). This
  function will return a two bit array structure allocated
  on the head with a two bit array of the requested size.
  The inputs are the number of elements in the two bit
  array and 0. However, if you wish to not make an array
  you can use twoBitStruct = makeTowBitArry(0, 1).

When finished the two bit array can be freed with
  the freeTwoBitAry function (fun-14 twoBitArrays.h). This
  takes a pointer to the twoBitAry structure to free,
  a 1 or 0 (free structure) to tell if to not free the
  structure, and a 1 or 0 (free array) to tell if to free
  the array stored in the structure. If you free the
  structure, make sure you set the pointer to null.

Their are several functions that are used to move around
  the two bit array, most of which only take a pointer to
  a two bit array structure or a pointer to the structure
  and a additional value.

To get the current two bit element on in the array use
  the getElmFromToBitUCAry() function
  (fun-01 twoBitArrays.h). To move forward one element use

You can get the index of the element you are on in the two
  bit array using the getTowBitAryIndex() function (fun-17
  twoBitArrays.h).

You can get the length of the two bit array by using the
  lenTwoBitAry() function (fun-16 twoBitArrays.c).
 
To change the current element on in a two bit array use the
  changeTwoBitElement() function (fun-08 in twoBitArrays.h)

You can do a shallow copy of the two bit array by using
  the cpTowBitPos() function (fun-13 twoBitArrays.h). This
  will allow you to work on the same two bit array with two
  separate pointers.
  
You can move forward in a two bit array using the
  twoBitAryMoveToNextElm() (fun-03 twoBitArrays.h) and
  twoBitAryMoveForXElm functions (fun-04 twoBitArrays.h).
  You can also use the moveXElmFromStart() function to move
  from the starting index instead of the current position.

You can move backwards in a two bit array using the
  twoBitAryMoveBackOnElm() (fun-05 twoBitArrays.h) and
  twoBitAryMoveBackXElm functions (fun-06 twoBitArrays.h).

Finally if you need to blank a limb to 0 you can use the
  blankLimb() function (fun-09 twoBitArrays.h).

Most of the movement functions have an alternative form
  that ends in BoundsCheck, which checks to see if the move
  resulted in an out of bounds error. I have not tested
  these, but they should work.

```
// An example of using two bit arrays

struct twoBitAry *twoBitST = makeTwoBitArry(10, 0);

// Add some values into the array
changeTwoBitElement(twoBitST, 1);
twoBitAryMoveToNextElm(twoBitSt);

changeTwoBitElement(twoBitST, 0);
twoBitAryMoveToNextElm(twoBitSt);

changeTwoBitElement(twoBitST, 3);
twoBitAryMoveToNextElm(twoBitSt);

// Move back to the start of the array
moveXElmFromStart(twoBitST, 0);

// Print out the first element
printf("%u\n", getElmFromToBitUCAry(twoBitST);

freeTwoBitAry(twoBitST, 0, 0);
twoBitST = 0;
return 0;
```
