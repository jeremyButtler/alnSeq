# Use:

This is here to give you an idea of the functions and data
  structures you will need to use to use this code in your
  own code.

This is not the best guide and does not cover everything,
  but I hope it will help.

## Note

See minExamples for a minimal C example. It currently is
  only for the Needleman without gap extension penalties,
  but should give you an idea on the flow of this code.

## Reading in sequences

Sequences are stored in the seqStruct, which is struct-01
  of seqStruct.h. seqStruct stores the sequence id,
  sequence, and q-score entries. It also stores the length
  of the id, sequence, and q-score and the size of each
  buffer storing the id, sequence, and q-score.

An seqStruct should be initialized with the initSeqST
  function (fun-07 in seqStruct.h) only for the first use.
  After initialization use the blankSeqST() function 
  (fun-08 seqStruct.h) to reset the structure to no
  sequence.

Use the freeSeqST function (fun-09 general/seqStruct.h)
  to clean  up a heap allocated seqStruct. freeSeqST takes
  a pointer to an seqStruct structure to free. Use
  blankSeqST (fun-08 general/seqStruct.h) to free an
  seqStruct that is on the stack (the sequence will be
  on the heap).

A sequence can be read in from a fastq file using the
  readFqSeq function (fun-02 in general/seqStruct.h) or
  fasta file using the readFaSeq function
  (fun-03 in general/seqStruct.h). Both functions take a
  pointer to the FILE with the sequence and a pointer to
  an seqStruct structure. The FILE pointer will be set to
  the next entry in the  fasta/fastq file.

If you compiled alnSetStruct.c without -DNOSEQCNVT, then
  you will also need to run seqToLookupIndex()
  (Fun-09 in general/alnSetStruct.h) to convert each base
  in a sequence to its index in the look up table. You can
  You can reconvert the sequence back to nucleotides with
  lookupIndexToSeq() (Fun-10 in general/alnSetStruct/h).
  Both lookupIndexToSeq and seqToLookUpIndex require a
  c-string (seqStruct.seqCStr), which is modified.
  Sll output bases will be in uppercase. These
  functions do nothing when -DNOSEQCNVT is used.

You can reverse complement a sequence using the
  reverseComplementSeq() function
  (fun-05 in general/seqStruct.h).

```
/*Here is an example.*/

FILE faFILE = fopen("sequence.fasta", "r");
struct seqStruct sequence;

if(faFILE == 0) return 0;

if(sequence == 0)
{
  fclose(faFILE);
  return 0;
}

initSeqST(&sequence);

switch(readFaSeq(faFILE, &sequence))
{ /*Switch: Check what kind of error*/
  case 0:
    fclose(faFILE);
    freeSeqSTStack(&sequence);
    sequence = 0;
    return 0;        /*EOF*/

  case 1:
    fclose(faFILE);
    break;
  case 2:
    fclose(faFILE);
    freeSeqSTStack(&sequence);
    sequence = 0;
    return 0;        /*invalid file*/

  case 64:
    fclose(faFILE);
    freeSeqSTStack(&sequence);
    sequence = 0;
    return 0;       /*memory error*/
} /*Switch: Check what kind of error*/

/*Only if aligning sequences*/
seqToLookupIndex(sequence.seqCStr);

/*Do alignment*/

lookupIndexToSeq(sequence.seqCStr);

/*Do something with the alignment*/

freeSeqSTStack(sequence);
```

## Alignment settings

Settings for the alignment settings are stored in the
  alnSet structure (struct-01 general/alnSetStruct.h).
  This structure holds the output line wrap size, which
  alignment to run, if you are doing a reference/query
  scan, the scoring matrix, the gap extension penalty, the
  gap starting penalty, the min score to keep an
  alternative alignment, and the order to prefer a single
  direction. The default values for these settings can all
  be modified in general/alnSeqDefaults.h.
  
You can initialize an alnSet structure using the
  initAlnSet() function
  (fun-11 in general/alnSetStruct.h). This will set all
  values in the alnSet structure to the default values in
  alnSeqDefaults.h.

For changing the scoring matrix you can use the
  setBasePairScore() function
  (fun-01 in general/alnSetStruct.h) to change the value
  for individual base pairs in a matrix or use the
  readInScoreFile() function
  (fun-07 general/alnSetStruct.h) to read in a scoring
  matrix. This scoring matrix is different than the normal
  format, so see scoring-matrix.txt for an example of my
  format.

For changing the match matrix (used in determing if had an
  match or a snp) you can use the setIfBpMatch() function
  (fun-02 in general/alnSetStruct.h) to change the value
  for individual base pairs in a matrix or use the
  readInMatchFile() function
  (fun-08 general/alnSetStruct.h) to read in a match
  matrix. This scoring matrix is different than the normal
  format, so see match-matrix.txt for an example of my
  format. Each match is treated as a 1 for is a match
  or 0 for is not a match.

Changing the non-scoring values for an alnSet structure
  requires you to manually change them. If the value ends
  in Bl, then it should only ever be a 1 or 0.

Finally you can free the alnSet structure that was put on
  the heap using the freeAlnSet() function
  (fun-04 in general/alnSetStruct.h). This takes in a
  pointer to an alnSet structure. Currently the alnSet
  structure does not allocate variables on the heap.
  However, there is a stack free (freeAlnSetStack)
  function (fun-03 general/alnSetStruct.h) which is
  present for future uses.

```
An example

struct alnSet settings;
FILE *matrixFILE = fopen("scoring-matrix.txt", "r");

if(matrixFILE == 0) return ERROR;

if(settings == 0)
{
  fclose(matrixFILE);
  return ERROR;
}

initAlnSet(&settings);
readInScoreFile(&alnSetST, matrixFILE);

/*Do something*/

fclose(matrixFILE);
freeAlnSetStack(&settings);

return 0;
```
  
## Doing an alignment

### General and default aligners

All alignment programs in alnSeq take three inputs,
  the query sequence as a seqStruct first, follwed by the
  reference sequence as an seqStruct, and then an
  alnSet structure with the alignment settings.

The return  values for each program can very. The
  Needleman and Waterman alignments return an
  alnMatrix structure if not doing alignments using or an
  two bit matrix or alnMatrixTwo bit structerfor
  alignments using a two bit matrix. The Hirschbeg
  returns an alnStruct, while the memory efficent
  waterman returns an alnMatrix structure with the
  index's for the alignment, but no matrix.

alnSeq can do either an Needleman Wunsch, Hirschberg, or
  an Waterman Smith alignment. The Needleman Wunsch
  alignment is done by the NeedlemanAln function
  (fun-01 in needleman/needleman.h). The Hirschberg
  is done by Hirschberg()
  (fun-02 hirschberg.h/hirschberg.h). The Waterman is done
  WatermanAln (fun-01 waterman/waterman.h). The memory
  efficent Waterman is the memWater() function
  (fun-01 memWater/memWater.h).

### Variations for each alignment algorithm

For the Waterman and memWater alignments you can also do
  a query reference scan; waterScan
  (fun-01 waterman/waterScan.h) and memWaterScan.h
  (fun-01 memWater/memWaterScan.h).

 For the Needleman you can also align without
  gap extension penalties (NeedleAlnNoGap(); fun-01 in
  needleman/needleNoGap.h), uses a two bit matrix
  (lower memory, but takes more time)
  (fun-01 needleman/needleTwoBit.h), and a Needleman that
  uses a two bit matrix and no gap extension penalty
  (fun-01 needleman/needleTwoBitNoGap.h)

For the Hirschberg you can also do an alignment without
  gap extension penalties with HirschbergNoGap (
  fun-02 in hirschberg/hirschbergNoGap.h).

For the Waterman you can do an alignment without gap
  extension penalties with WatermanAlnNoGap()
  (fun-01 in waterman/watermanNoGap.h), with two bit
  arrays with waterTwoBit() (fun-01 in
  waterman/waterTwoBit.h), and with two bit arrays and
  without gap extension penalties with WaterTwoBitNoGap()
  (fun-01 in waterman/waterTwoBitNoGap.h).

For the query reference Wateman scans you can also not use
  gap extension penalties with WaterScanNoGap (fun-01
  waterman/waterScanNoGap.h), using two bit arrays with
  WaterScanTwoBit (fun-01 waterman/waterScanTwoBit.h), and
  with two bit arrays and no gap extension penalites
  with waterScanTwoBitNoGap() (fun-01
  waterman/waterScanTwoBitNoGap.h).

All alignment functions take an pointer to an seqStruct
  structure with the reference sequence, an pointer to an
  seqStruct structure with the query sequence, and an
  pointer to a alignment settings structure.

The Needleman Wunsch, Waterman Smith, and query/reference
  scan waterman alignments will return an alnMatrixStruct
  structure (struct-01 alnMatrixStruct.h). The
  alnMatrixStruct has the directional matrix,
  best score (Waterman)/ending score (Needleman), and the
  extra scores kept for the reference/query scan
  (if was called for).

For the memory efficent Waterman you can do an alignment
  that does not use gap extension penalties (only gap
  open) with memWaterNoGap (fun-01
  memWater/memWaterNoGap.h). You can also do an alignment
  with gap extension penalties for the memWater with
  memWaterScanNoGap() (fun-01
  memWater/memWaterScanNoGap.h).

## What to do after the alignment

The Needleman, Waterman, and memory efficient Waterman all
  the end of the alignment
  (matrix->bestEndIndexUL). They also return the best
  score (matrix->bestScoreL). However, the Needleman and
  Waterman also return a directional matrix
  (matrix->dirMatrix) which can be converted to an
  alnStruct (struct-01 general/alnStruct.h) by using
  the ending index with dirMatrixToAln()
  (fun-05 general/alnStruct.h) and the ending index
  position or twoBitDirMatrixToAln()
  (fun-06 general/alnStruct.h) for matrix's that used two
  bit arrays. This alnStruct can then be used to get an
  alignment with alnSTToSeq (fun-04 general/alnStruct.h)
  or prnted with printAln (fun-15 general/alnStruct.h).

The memory efficent Waterman will also return a starting
  index (matrix->bestStartIndexUL). The starting and
  ending index's can be converted to the rerence and query
  position using indexToCoord() (fun-09 in
  general/alnMatrixStruct.h). This is a macro that takes
  in the reference length (matrix->lenRefUL), index,
  variable to hold the reference coordinate (do not use
  pointers), and the variable to hold the query coordinate
  (again do not use pointers).

The query reference scans also modfiy the alnMatrix or
  alnTwoBitMatrix to hold the scores for each kept base
  (matrix->scoreAryL), the starting index of each kept
  score (matrix->startIndexAryUL), the ending index of
  each kept score (matrix->endIndexAryUL), and a variable
  telling the number of maximum kept scores
  (matrix->lenArraysUL).

Once you are finished with the matrix you can free it with
  the freeAlnMatrix() function (fun-05 in
  general/alnMatrixStruct.h) or for matrixs with two bit
  arrays use freeAlnMatrixTwoBit() (fun-06 in
  gernal/alnMatrixStruct.h).

```
struct seqStruct querySeq;
struct seqStruct refSeq;
struct alnSet settings;

struct alnMatrix *matrix = 0;
struct alnStruct *alignmentST = 0;

FILE *faFILE = 0;

initAlnSet(&settings);
initSeqST(&querySeq);
initSeqST(&refSeq);

faFILE = fopen("query.fasta", "r");
if(faFILE == 0) return 0;
fclose(faFILE);

switch(readFaSeq(faFILE, &querySeq))
{
/*Check to see if read in sequence & handle errors here*/
}

faFILE = fopen("ref.fasta", "r");

if(faFILE == 0)
{
/*Handle file error here*/
}

fclose(faFILE);

switch(readFaSeq(faFILE, &refSeq))
{
/*Check to see if read in sequence & handle errors here*/
}

queryST.endAlnUI = queryST.lenSeqUI;
queryST.offsetUI = 0;

refST.endAlnUI = refST.lenSeqUI;
refST.offsetUI = 0;

matrix = WatermanAln(&querySeq, &refSeq, &settings);
alignmentST =
    dirMatrixToAln(
       &refST,
       &queryST,
       matrix->bestEndIndexUL,
       &settings,
       matrix
);

/* or matrix = NeedlemanAln(&querySeq,&refSeq,&settings);
` or alignmentST = Hirschberg(&querySeq,&refSeq,&settings);
*/

/* Memory efficent waterman
`  bestScoreST = memWaterAln(&queryST, &refST, &settings);
`  if(bestScoreST == 0) 
`  {
`   // Handle Errors
`  }
`
`  indexToCoord(
`     matrix->lenRefUL, 
`     matrix->bestStartIndexUL,
`     refST.offsetUL,
`     queryST.offsetUL
`  );
`  indexToCoord(
`     matrix->lenRefUL, 
`     matrix->bestStartIndexUL,
`     refST.endAlnUL,
`     queryST.endAlnUL
`  );
`
`  alignmentST = Hirschberg(&refST, &queryST, &settings);
*/

/*Do something with alignment here*/

if(bestScoreST != 0) freeScoresST(bestScoreST, 0);
if(matrix != 0) freeAlnMatrixST(matrix);
if(alignmentST != 0) freeAlnST(alignmentST, 1);
freeSeqST(&refSeq, 0);
freeSeqST(&qrySeq, 0);
freeAlnSet(settings, 0);
```
  
## Printing alignments

### Printing alternative alignments from a scan

Alternative alignments can be printed with the
  printAltWaterAlns function (fun-01 general/genScan.h).
  This function takes the alnMatrix or alnMatrixTowBit
  structure returned by WatermanScan and memWaterScan,
  a min score to keep an alternative alignment, and an
  output file.

### Printing primary alignments

To print out the alignment you first have to create an
  alnStruct. This is already done for the Hirschberg 
  alignment, but not for the Needleman or Waterman
  alignments. 

Converting the directional matrix to an alignment array 
  is done with the dirMatrixToAln) or towBitDirMatrixToAln
  (fun-05 and fun-06 in general/alnStruct.h). These
  functions take in a pointer to an seqStruct structure
  with the reference sequence, an pointer to an seqStruct
  structure with the query  sequence, the index of the end
  of the alignment (matrix->bestEndIndexUL or
  matrix->endIndexAryUL[position], and the directional
  matrix returned from the Needleman or Waterman
  alignment.

After building the alignment array you can then free the
  alnMatrixStruct structure. Just remember to save the best
  score before freeing the alnMatrix structure.

You can print out the alignment stored in the alnStruct
  using the printAln function (fun-15 in
  general/alnStruct.h).

Just make sure to free the alignment array structure
  (alnStruct) with freeAlnST(alnStruct);.

```
char *alignedQuery = 0;
char *alignedRef = 0;
long bestScoreL = 0;
struct alnStruct *alignment = 0;

/*Other variables, such as the output file*/

seqToLookupIndex(refSeqST.seqCStr);
seqToLookupIndex(qrySeqST.seqCStr);

matrix = WatermanAln(refSeqST, querySeqST, settings);

/*Get the alignment array*/
alignment =
  dirMatrixToAlnST(
    refSeq,
    querySeq,
    matrix->bestEndIndexUL,
    matrix
); /*Build the alignment array*/

bestScoreL = matrix->bestScoreL;

freeAlnMatrix(matrix);
matrix=0;

if(alignment == 0)
{
  /*Handle memory errors*/
}

/*Printing out the structure*/

lookupIndexToSeq(refSeqST.seqCStr);
lookupIndexToSeq(querySeqST.seqCStr);

printAln(
  outFILE,              /*File to output to*/
  outFileStr,           /*File name*/
  refSeqST,
  querySeqST,
  alignment 
  bestScoreL,
  settings,             /*unsigned short from alnSetStruct*/
  alignment             /*alnStruct*/
  scoreMatrixFileNamea  /*Name of user input score matrix*/
    /*Use defMatrixNameStr if using the default matrix*/
);

fclose(outFILE);

freeSeqST(refSeqST);
freeSeqST(querySeqST);

freeAlnST(alignment, 1);
freeAlnSet(settings, 0);

return 0;
```

## Two bit arrays

All twoBit functions are static inlined (only a
  twoBitArrays.h file).

The two bit array is an array of uint8_t's that each hold
  four two bit elements. This array is stored in a
  twoBitAry structure (struct-01 twoBitArrays.h). The
  twoBitAry structure has a pointer to the first unit8_t
  (limb) in the array (firstLimbUCPtr), a pointer to the 
  limb currently one (limbOnUCPtr), a counter telling which
  element we are on in the limb, and the length of the
  array.

You can create a twoBitAry structure using the makeTwoBit
  function (fun-12 in twoBitArrys.h). This function will
  return a two bit array structure allocated on the head
  with a two bit array of the requested size. The inputs
  are the number of elements in the two bit array and 0.
  However, if you wish to not make an array you can use
  twoBitStruct = makeTowBit(0, 1).

When finished the two bit array can be freed with
  the freeTwoBit function (fun-14 twoBitArrays.h). This
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
  the getTwoBitElm() function
  (fun-01 twoBitArrays.h). To move forward one element use

You can get the index of the element you are on in the two
  bit array using the twoBitIndex() function (fun-17
  twoBitArrays.h).

You can get the length of the two bit array by using the
  twoBitGetLen() function (fun-16 twoBitArrays.c).
 
To change the current element on in a two bit array use the
  changeTwoBitElm() function (fun-08 in twoBitArrays.h)

You can do a shallow copy of the two bit array by using
  the cpTwoBitPos() function (fun-13 twoBitArrays.h). This
  will allow you to work on the same two bit array with two
  separate pointers.
  
You can move forward in a two bit array using the
  twoBitMvToNextElm() (fun-03 twoBitArrays.h) and
  twoBitMvForXElm functions (fun-04 twoBitArrays.h).
  You can also use the twoBitMvXElmFromStart() function
  (fun-15 twoBitArrays.h) to move from the starting index
  instead of the current position.

You can move backwards in a two bit array using the
  twoBitMvBackOnElm() (fun-05 twoBitArrays.h) and
  twoBitMvBackXElm functions (fun-06 twoBitArrays.h).

Finally if you need to blank a limb to 0 you can use the
  blankTwoBitLimb() function (fun-09 twoBitArrays.h).

```
/*An example of using two bit arrays*/

struct twoBitAry *twoBitST = makeTwoBit(10, 0);

/*Add some values into the array*/
changeTwoBitElm(twoBitST, 1);
twoBitMvToNextElm(twoBitSt);

changeTwoBitElm(twoBitST, 0);
twoBitMvToNextElm(twoBitSt);

changeTwoBitElm(twoBitST, 3);
twoBitMvToNextElm(twoBitSt);

/*Move back to the start of the array*/
twoBitMvXElmFromStart(twoBitST, 0);

/*Print out the first element*/
printf("%u\n", getETwoBitElm(twoBitST);

freeTwoBit(twoBitST, 0, 0);
twoBitST = 0;
return 0;
```