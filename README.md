# Use

AlnSeq uses a Smith Waterman and Needleman Wunsch alignment
  that runs with a memory usage of (O(n \* m / 4)).
  However, it is less speed efficient than a traditional 
  Smith Waterman and Needleman Wunsch alignment. It can
  not do matrix scanning, but will allow for multiple local
  alignments to be output from a single alignment.

One thing I should add. I am not well read on the alignment
  literature, so I have no idea if some one has already
  done this. My impression is no, because I do not see
  anything like this done in emboss, mentioned on
  Wikipedia, or mentioned in the two or three reviews I
  have read on the Smith Waterman algorithm.

# This code is not running yet

I am currently debugging this code and will hopefully be
  done in a few weeks to months. For an alternative,
  working version of this code please see the alingSeq
  folder at
  [https://github.com/jeremybuttler/find--Co-infections](https://github.com/jeremybuttler/find--Co-infections).
  The only major difference between the alternative code
  and this code will be that this code will allow the
  option of saving more then one alignment for the Smith
  Waterman alignment.

Here is my current status

- The Needleman Wunsch alignment has be tested and is now
  working on OpenBSD (the OS I mainly use), but it is only
  working on the Debain Linux machine I am testing with
  non-stdout output. Use a ">", "|", or "-out file.aln" to
  no error out.
- The Smith Waterman alignment is working, but printing
  out needs a redirect like the Needleman Wunsch on linux.
- I am not sure if multi-output Smith Water alignment is
  working. I will start debugging this next.

# Likely better alternatives already on github

- https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library

# Building and running alnSeq

## How to build alnSeq

```
# How to install this program

# Install at /usr/local/bin (need root privilage)
make
sudo make install

# Use a different install path
make
make install PREFIX=/path/to/install

# Manual install
make
mv alnSeq /path/to/install
chmod a+x /path/to/install/alnSeq
```

## How to run alnSeq

```
# For the help message
alnSeq -h

# For a global alignment (Needleman Wunsch)
alnSeq -query query.fasta -ref ref.fasta > alignment.aln

# For a single local alignment (Waterman Smith) (not working)
alnSeq -use-water -query query.fasta -ref ref.fasta > out.aln

# For a multi local alignment (Waterman Smith) (not working)
alnSeq -multi-base-water -query query.fasta -ref ref.fasta -out prefix 

```

# Explaining alnSeq

## What is alnSeq?

AlnSeq is an sequence alignment program that uses either 
  a Needleman Wunsch or Smith Waterman alignment algorithm.
  AlnSeq is unique from a traditional Needleman Wunsch or
  Smith Waterman in how it handles the scoring matrix and
  the direction matrix.

One thing I do want to point out is that there are very
  memory efficient algorithms for optimal global alignments
  (Needleman Wunsch), such as the Hirschberg alignment.
  However, I am unaware of any solutions for optimal local
  alignments (Smith Waterman). I am not very well read in
  the literature, so it is possible I missed something.

## The direction matrix

AlnSeq stores each direction in two bits, which are packed
  into an 8 bit integer. This reduces the directional
  matrix size by 4, but comes at the cost of slower speed.
  This also means that only one direction is stored,
  instead of all possible alternatives. However, alnSeq
  could be modified to keep all directions in a 3 or 4 bit
  array was used at the cost of using twice the memory. You
  can modify the twoBitArrays.c functions to make this
  change.

If you do decide to store all directions , then I would
  recommend using 3 bit directions packed into a 16 bit or
  32 bit integer. It will be harder and will loose one bit
  per 16bit integer, but will increase memory usage by 60%
  instead of 100%.

## The scoring matrix

AlnSeq also reduces the scoring matrix down to two rows,
  which are the previous row of scores and the current row
  of scores. This effectively removes the scoring matrix,
  but also removes the ability to scan the scoring matrix
  for other high scores, which is sometimes done for a
  traditional Waterman Smith alignment. To combat this
  alnSeq will record the best score (or lower right conner
  for a Needleman Wunsch) and the best score for each base
  in the reference and each base in the query sequence.

To reduce the tendency to trace the path of the best
  alignment each new best score for a bases is check to
  see if it extends the alignment the previous base was on.
  If it extends the alignment, then the previous bases best
  score is set to an previous score for the previous base
  and the current base is set to the new score. However,
  this means that alnSeq must keep track of both the
  current best scores for each base and one past best
  scores for each base.

Any bases that had a paths (alignments) that had a high
  enough score are printed out. The total number of printed
  alignments can be thought of in terms of n (number of 
  bases in the reference) and m (number of bases in the
  query). The maximum of number of alignments is n
  alignments for the reference and m alignments for the
  query. This could take up a lot of space, so it would
  be better to integrate your program into this step. You
  can use the printAlnWaterAlns (Fun-03) in water.c, which
  does the printing as an example. 

From this you could get an idea of the complete alignment
  and may even detect duplications that are in only one 
  sequence. However, detecting duplications that are in
  both sequences will likely not be possible.

If your willing to modify the code you could make an
  alternative system that stores every score above a
  certain threshold or multiple scores per base. For the
  every score above a certain threshold modification you
  would need to remove the old scores arrays, which when
  done would allow you to store up to 25% of the matrix
  before having similar memory usage to a traditional
  Smith Waterman. If you did not care about scores you
  could increase this to 50% of the matrix.

## Memory Cost

The rough planned memory cost of alnSeq is

- Cost of the directional matrix (n\*m)/4
- Cost to store the reference and query sequence (n + m)
- Cost of storing the best index 8 \* (n + m) \* 2
  - The index must be stored as unsigned longs
  - \*2 is because the previous best score needs to be
    stored
  - This is for each query and reference base, but not the
    entire matrix
- Cost of storing every best score 8 \* (n + m) * 2
  - Scores are stored as longs
  - \*2 is because the previous best score needs to be
    stored
  - This is for each query and reference base, but not the
    entire matrix
- Cost of converting matrix to an alignment 3 \* (n + m)
- Total rough cost (n\*m)/4 + 36 \* (n + m)
- Cost of keeping one more score per base 32 \* (n + m)
  
The memory cost for two 200,000 base pairs sequences would
  be ([(200,000 \* 200,000) / 4] + 36 \* 200,000) bytes =
  10 Gb + 7.2 Mb < 11Gb.

A traditional Smith Waterman alignment would take 4 bytes
  per cell of the scoring matrix and at best could take
  1/4 of a byte for the directional matrix, but often takes
  1 byte. This means that for aligning two 200,000 base
  pair sequences the traditional Smith Waterman uses
  4 to 5 bytes \* (200,000 \* 200,000) bytes = 160Gb-200Gb.
  This is large increase in memory usage, but also allows
  the detection of all duplications with a matrix scan.

This shows that though memory may be cheap, it is not n^2
  cheap.

## Final notes

AlnSeq does not use decimals, so if you want decimals for
  the gap extension penalty, so you will have to multiply
  all scores by 10 to 1000.

My current vision for this code is more to be an example or
  demonstration of an idea, rather then to have all
  features added in. This means that I will be done with
  alnSeq when I am done debugging and that I do not plan to
  take alnSeq any farther right now. This means some of the
  additions I mentioned, multi-threading, or even gpu
  support is up to you.

Note: This program is more set up for noisy reads. This
  means that the gap opening has a lighter penalty than the
  gap extension penalty, which means that this program
  favors single gaps instead of longer gaps. For a
  reference or consensus alignment you will want to set
  the gap open penalty higher and lower the gap extension
  penalty.

For multi-threading my thought was to do rows and have each
  thread launch a thread after it finishes the first cell
  in the row. I then planned to keep an additional row to
  have each thread mark their row number for each score
  they complete. The idea is that a thread only advances to
  the next cell in the matrix when the required previous
  scores have been found and marked from the previous row.
  To conserve memory I would have the threads overwrite the
  previous scores in the scoring matrix for each cell they
  complete. You would need a lock when recording the best
  score(s).

This would take some work, but you could do a full scan of
  the directional matrix if you recalculated all the
  scores. However, this will likely take more time than a
  traditional Smith Waterman. This should take less time
  than the original calculations though because, you only
  need to find one score per cell and already know the
  chosen direction for each cell.

1. Check direction
  - If match/snp compare bases to get score
  - Otherwise check previous direction to see if applying
    an indel starting or indel extension penalty
2. Get score of cell by adding match/snp or the indel
   penalty to the previous this cell pointed to (the cell
   in its path)
3. If the score is within range, then skip rest of row and
   complete the alignment path
   - use a separate scoring array for this
   - Ignore all cells not in the path of this alignment
4. Once the path is complete, print, check, or doing what
   you want with the alignment
   - Path is complete once at a stop or hit a cell with
     a lower then minimum score
   - This may require you to skip some alternative paths
5. Mark all cells in the traveled path as a stops, so this
   path is not traveled again
   - As you come up, travel down all missed valid
     alternative paths and repeat 3 to 5 till no
     alternative paths exist.
6. Return back to start of path to complete the row/skipped
   cells, repeating steps 1 to 6 until all cells have been
   checked.

# Thanks

- To my dad for being willing to listen and give advice.
- To the people (never met) who coded baba
  [https://baba.sourceforge.net/](https://baba.sourceforge.net/).
  It gave me a great visual on how the Needleman Wunsch
  algorithm worked. Their were many other sources, but this
  was the one that was the most useful to me.
