# Use

This log records how alnSeq has changed between versions.
  Each version is labeled with the version number, year,
  month, and then day. For example version 1.20230615 would
  be version 1, from June 15th, 2023.

# Issues to be solved

- multi sequence output for Waterman Smith alignment
- Printing to stdout on Debain, which will likely be
  similar for all linux's


# Issues that I am not sure if I will work on

- Full matrix scan for Waterman Smith
- Multithreading & gpu support with opengl

# Log

## Version 1.20230615

- Smith Waterman single answer alignment is now working
- Fixed some issues with printing out alignments.
  - Converting matrix index's to base positions was off
  - Alignment printing function was printing extra bases
    or two few bases

## Pre Version 1.20230615

- Logic set up and Needleman Wunsch alignment working
