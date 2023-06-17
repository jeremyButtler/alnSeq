# Use

This log records how alnSeq has changed between versions.
  Each version is labeled with the version number, year,
  month, and then day. For example version 1.20230615 would
  be version 1, from June 15th, 2023.

# Issues to be solved

- multi sequence output for Waterman Smith alignment
  - Currently segfaults
- Matrix scan for Waterman Smith
  - Reference starting base is going out of bounds
- Printing to stdout on Debain, which will likely be
  similar for all linux's
  - This is likely an array out of bounds error some were.
    I suspect this is in my convert directional matrix
    function to an alignment array.

# Ideas that would be cool, but not worth working on

These ideas are a future vision that is not worth the
  time unless this algorithm becomes needed by others.

- Multithreading & gpu support with opengl
- Converting Smith Waterman matrix to a jagged cigar
  matrix.
  - This would use a three bit matrix, with the third
    bit marking if the next 16 bit int is a number.
  - This would be slower, but also be more compressed.

# Log

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
