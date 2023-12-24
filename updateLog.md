# Use

This log records how alnSeq has changed between versions.
  Each version is labeled with the version number, year,
  month, and then day. For example version 1.20230615 would
  be version 1, from June 15th, 2023.

# TODO or to fix list

- Vector support for memory efficient Waterman and
  Hirschberg.
- Fix mem water giving different anwser for reverseing
  the reference and query. At least for the small
  test using the query as the reference ignores the
  the starting insertions.
- Python wrapper returns a full alignment when doing
  paritial alignments. This is also true for the
  C code when specifing a specfic range to align
- Get the filters for query reference scans working

# Ideas that would be cool, but not worth working on

These ideas are a future visions that are not worth the
  time unless this algorithm becomes useful to others.

- Multithreading & gpu support with openCL
- Get a matrix scan function that recalculates scores on
  fly, so that user can search a completed direction
  matrix.
- Add a filter for alternative alignments.

# Log

## 20231224

- Fixed some small errors that were resulting in random
  indels when alining larger genomes
- Removed all .c files, except the main file. There are
  some functions that are set with "static" and so will
  raise the -unused-function flag in gcc. Use the
  -Wno-unused-function flag when compiling.
- Converted loops for alignments into for loops. This gave
  a minor speed boost for all functions.
- Fixed query reference scans, and changed the requirement
  to keep a score to be a snp. I also make sure each score
  is only recoreded once. For the frist half of the
  reference it is the reference scores, then thw query.
  For the last half  it is first recoreded under the query
  scores, then the reference.
  - This change has made the scans a bit slower.
- I am now recording index's instead of query and
  reference starting positions. This means the maximum
  alignment size for the memWaters is the maximum size
  of an unsined int (~ 4 billion) - 1.
  - This has speed up the memory efficent waterman
    alignments.
  - It also would speed up the query reference scans,
    except I changed how things are recoreded.
- I have set up a separate function for normal alignments,
  two bit alignments (-two-bit) (waterman/needleman only),
  and for  alignments without gap extension penalties
  (-no-gapExtend).
  - For python this has been added in the wrapper
    functions as noGapBool (no gap extension) and
    twoBitBool (use to bit arrays).
- Added in the ability to set what counts as a match.
  See the match-matrix.txt file for an example.
- alnSeq is now c89 complient
- Updated usingThisCode.md to reflect the new changes.
  Not everything will be in there, but it should be
  a start.

## 20231029

- Changed complile settings for alnSeq to have -DBYTEMATRIX
  and -DINSDELSNP. These can be disabled by overwriting
  the CFLAGS or CFLAGS="-DBLANK".
  - This is the old mid compile option.
- Set up alnSeq to be a library in Python (pythonPkg).
  - Supported: Hirschberg, Needleman Wunsch, and
    Smith Waterman
- Fixed a minor bug in the Hirschberg conversion to an
  alnStruct function. In rare cases it caused an infinite
  loop.
- Add in python library compile option.

## 20231022

- Changed the default settings for alnSeq to favor
  insertions, deletions, and then snps/matches. This was
  to deal with an error I caught in some of my alignments.
- Set the scoreing step in the Hirschberg to only recored
  if direction was a gap or snp (instead of ins, del, or
  snp). This speeds up the Hirschberg at no cost to
  accuracy.
  - This change broke the twobit compiled Hirschberg
    (-DHIRSCHWOBIT). My best guess it this is an error in
    my if defines.
- Changed how I found the max direction in the Needleman,
  Waterman, and memory waterman programs. This reduced the
  number of operations and provided a small speed increase.
  - Old way was a maximize function.
  - This change broke the query scans and twobit memory
    efficent Smith Waterman (-DTWOBITMSW). I suspect that
    this may be a copy/paste error, but I could be wrong.

## 20230908

- Fixed an issue in the memory efficent Watermans were the
  alternative alignments were not printing out the starting
  alignments.
  - This is still a bit broken for the non-memory efficent
    Waterman. This was due to an issue with the Hirschberg
    adding bases at the end for stich (my
    artic-benchmarking repository). I have been unable to
    replicate this error, so it may have been an issue
    with stich.

## 20230831

- Fixed insertion error for memory efficient Waterman.

## 20230830

- Added in the memory efficient Waterman.
- Changed -query-ref-scan-water to -query-ref-scan.
  - The memory efficient Waterman is now the default option
    when -query-ref-scan is used. This can be changed by
    using -use-water.
- Fixed a minor error in my if not defines in waterman.c,
  which resulted in a byte matrix being used in the
  alternative alignment function with -DNOGAPOPEN.

## 20230827

- Alternative alignment printing now prints out the
  score, starting bases of alignments, and ending bases of
  alignments (Used to be cigars).
  - Still need to add a filter.
- The query-reference scan has been simplified.
- Matrix scan has been removed. It was not that great in
  the first place.
- Issue with Waterman not printing out full length
  alignments has been fixed.
  - It was do to using a `&&` when I should have been
    using `||` in the while loop in printAln
    (alnStruct.c Fun-03 Sec-05 Sub-01).
    ```
    *refAlnStr != defEndAlnFlag ||
    *qryAlnStr != defEndAlnFlag
    ```
- Fixed an issue in -format-emboss, were it would print
  out the wrong alignment length. It now prints out the
  matches + insertions + deletions instead of the longest
  sequence.
- Adding in several compiler time flags
  - -DNOGAPOPEN: Disables gap opening penalties
    - Default gap penatly is set to -6 for this option
    - Removed -gapopen option.
  - -DBYTEMATRIX: Use a character matrix instead of an
    two bit matrix for the Needleman and Waterman
    alignments.
  - Flags to make alnSeq use only one direction
    - Removes the direction selection options.
    - -DSNPINSDEL
    - -DSNPDELINS
    - -DINSSNPDEL
    - -DINSDELSNP
    - -DDELSNPINS
    - -DDELINSSNP
- Added in a `-flag` and `-flag-only` options to print out
  the flags alnSeq was complied with. `-flag` will also
  print out a description of each flag.
- Shortened two bit array function names
  - TwoBitArry -> TwoBit
    - `sed 's/TwoBitArr*y/TwoBit/g'`
  - TwoBitAry  -> TwoBit
    - `sed 's/TwoBitArr*y/TwoBit/g'`
  - twoBitAryMove -> twoBitAryMv
    - `sed 's/twoBitAryMove/twoBitAryMv'`
  - getTwoBitAryIndes -> twoBitGetIndex
- Added to twoBit to two bit function names missing two bit
  - moveToNextLimb    -> mvToNextTwoBitLimb
  - moveToLastLimb    -> twoBitMvToLastLimb
  - moveXElmFromStart -> twoBitMvXElmFromStart
  - lenTwoBit         -> twoBitGetLen
- Updated using this code guide for the new printing and
  alnStruct changes.

## Version 20230811

- Switched from converting my bases to a look up table on
  checks to the read in sequence.
  - (disable with `make CFLAGS="-DNOSEQCNVT"`)
- Improved speed by inlining several functions and adding
  some branchless operations at key points
  (selecting max score, and zeroing waterman values). Other
  places not worth it.
  - All functions in twoBitArrays.c have been inlined and
    are now in twoBitArrays.h.
  - Scoring functions are macros in generalAlnFun.h
- Hirschberg is now compiled with 1 byte arrays by default
  - This is faster, but use a small amount of extra memory.
  - enable twobit version with
    `make CFLAGS=`-DHIRSCHTWOBIT`.
- Made Hirschberg thread safe. This will increased memory
  usage by a very slight amount
  - ~ bytes = 1/4 \* reference length for two bit.
  - ~ bytes = \* reference length for byte arrays
  - This was added in at this point, but I forgot to 
    mention the updated. That is why it is added in at
    a later date.


## Version 20230804

- Fixed base getting removed on printing. It turned out to
  be error in my fasta reader (assigning an index 0
  variable an index 1 value). So, it was not a printing
  error.
  - line 510 in seqStruct.c `spareBuffUL = resBuffUL;` was
    changed to `spareBuffUL = resBuffUL - 1;`.

## Version 20230803

- The Hirschberg now prints using the normal functions.
  - I corrected an error were I was sending in the query as
    the reference to the Hirschberg.
- alnStruct has been changed. It now stores two alignment
  arrays. One for the reference and one for the query.
  - Each array tells if a specific reference base or query
    base mapped to a gap, snp, match, or was soft masked.
  - Also added in are insertion, deletion, snp, and match
    counts.
  - Also includes the first match/snp and last base in the
    alignment
  - The alignment arrays now include matches.
- The Smith Waterman is now printing correctly when a line
  wrap is applied.
- Changed printAln (alnStruct.c/h, function 03, used to be
  04) to take in more input. This allowed me to print
  out more to the header.
  - You can now print the output file in fasta, clustal,
    or an emboss like format. The old expanded cigar
    format is still the default option.
  - I am hoping the change in printing method got rid of a
    rare error I found a few days ago, were one or two
    bases would be replaced with an deletion. I have not
    confirmed if my hope is correct.
- Added in an option `-print-alinged` to only print out the
  aligned bases in the alignment.
- Added in option `-print-position` to include starting and
  ending base numbers at each line.
  - Always set when output is in the EMBOSS like format.
- Updated the help message. I also made a printHelpMesg
  function, which adds in the default values in 
  alnSeqDefualts.h/c to the help message and prints it out.
  This should avoid the help message getting out of date
  in the future.

## Version 20230726

- Hirschberg is working, but the output methods are
  limited.
- Gap opening penalty has been changed to -10. I found this
  worked better when using alnSeq for my projects.

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
