# Use

This log records how alnSeq has changed between versions.
  Each version is labeled with the version number, year,
  month, and then day. For example version 1.20230615 would
  be version 1, from June 15th, 2023.

# Ideas that would be cool, but not worth working on

These ideas are a future vision that is not worth the
  time unless this algorithm becomes useful to others.

- Multithreading & gpu support with openCL
- Converting Smith Waterman matrix to a jagged cigar
  matrix.
  - This would use a three bit matrix, with the third
    bit marking if the next 16 bit int is a number.
  - This would be slower, but also use a more compressed
    alignment matrix.
- Get a matrix scan function that recalculates scores on
  fly, so that user can search a completed direction
  matrix.
  - Also allow a function to be input, so user can decide
    if printing or not.

# Log

## Version 20230709

- using "-line-wrap 0" will now disable line wrapping
- Worked out some bugs in the Hirschberg alignment
  - Hirschberg still not working

## Version 1.20230630

- Changed sequence lengths and offsets to be longs instead
  of integers
- Put down logic for Hirschberg, which I am currently
  debugging
- Changed the default gap extension and gap opening 
  penalties to be more suited for non-noisy reads
  - Gap opening penalty is now -4 (used to be -1)
  - Gap extension penalty is now -1 (used to be -4)
- Changed
  
## Version 1.20230620

- Fixed an issue were two very different large alignments
  (Tested was: testing/largeTest-query-a.fasta and
  TickBornEncephalitis-reference.fasta) would segfault.
  - This error was due to my inputting unsigned longs
    instead int32_t for twoBitAryMoveForXElm.
  - twoBitAryMoveForXElm now uses unsigned longs instead of
    int32_t for its index

## Version 1.20230619

Other than bug fixes that pop up, this will be my final
  code update until I see that this code is useful to
  others. The only exception is if I find a use for 
  pressing this code further.

- Fixed the printing to stdout on Debain error.
  - This was due to alnST not being allocated properly
    in the cnvtDirMatrixToAlnAry function with calloc. This
    resulted in only a partial alnStruct structure being
    allocated.
  - Also had problem were I freed the alignment matrix to
    soon and so lost the best score. This has been fixed
    by storing the best score in a temporary variable.
- Multi sequence output that prints a best score for each
  query and reference for Waterman Smith alignment no
  longer segfaults
  - Was problem with using double pointers to make an array
    of structures.
  - Setting has been changed to "-query-ref-scan-water".
- Fixed matrix scan for Waterman Smith having reference
  starting base is going out of bounds of unsigned integer.
  - This was caused by my Loop always going negative for a
    counter, so I switched to finding starting base by
    using the ending index (end of the traced alignment
    path or start of alignment) of the two-bit Array to
    find the starting base.
- The reported base range for the alignment is now in index
  0 (first base is position 0) instead of index 1 (first
  base is position 1).
- The min score for matrix scanning and query/reference
  scanning (best score for each reference/query base) has
  been changed to 1000 (at least 200 bases).
  - Score increased for each non-anonymous match is 5.
- Set U and T to be separate scores and updated the default
  EDNA scoring matrix to have U and T be the same.
  - This was set because I built alnSeq to align DNA and
    RNA sequences, which have T and U being the same. This
    would have had no effect on protein alignments, which
    do not use U, but would have had an effect on word
    alignments.
  - AlnSeq ignores case for all alignments.
  - The scoring matrix alnSeq can be set up to have scores
    for the characters a-z.
    
## Version 1.20230616

- The Smith Waterman multi base alignment logic has been
  modified to only record a best score if the alignment
  does not extend to the next row or the current score
  is better then the next score in the alignments path.
  - Currently this segfaults
- Added in a matrix scan feature.
  - This currently is off on the reference starting base

## Version 1.20230615

- Smith Waterman single answer alignment is now working
- Fixed some issues with printing out alignments.
  - Converting matrix index's to base positions was off
  - Alignment printing function was printing extra bases
    or two few bases

## Pre Version 1.20230615

- Logic set up and Needleman Wunsch alignment working
