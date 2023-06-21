# Use

AlnSeq uses a Smith Waterman and Needleman Wunsch alignment
  that runs with a memory usage of (O(n \* m / 4)).
  However, it is less speed efficient than a traditional 
  Smith Waterman and Needleman Wunsch alignment. It can
  output multiple local alignments from a single alignment.

There are faster and less memory hungry Waterman Smith
  implementations than alnSeq. One example is the stripped
  Waterman Smith alignment, which I think reduces both
  scoring and direction matrix to just a few rows. See 
  [https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
  for an very fast striped Waterman Smith aligner.

Other than doing some bechmarking and future debugging, I
  am finished with this code. If this code becomes usefull
  to many others I may try to improve this code further,
  but for now I am done.

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

# For printing out a best score for each base
alnSeq -query-ref-scan-water -query query.fasta -ref ref.fasta -out prefix 

# For a matrix scan
alnSeq -matrix-scan-water -query query.fasta -ref ref.fasta -out prefix 

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
  The striped Waterman Smith alignment is a very memory
  efficient algorithm for optimal global alignments.

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
- Cost of storing the best index 8 \* (n + m)
  - The index must be stored as unsigned longs
  - \*2 is because the previous best score needs to be
    stored
  - This is for each query and reference base, but not the
    entire matrix
- Cost of storing every best score 8 \* (n + m)
  - Scores are stored as longs
  - This is for each query and reference base, but not the
    entire matrix
- Cost of converting matrix to an alignment 3 \* (n + m)
- Total rough cost (n\*m)/4 + 19 \* (n + m)
- Cost of keeping one more score per base 16 \* (n + m)
  
The memory cost for two 200,000 base pairs sequences would
  be ([(200,000 \* 200,000) / 4] + 19 \* 200,000) bytes >
  10 Gb, but < 11Gb.

This shows that alnSeq has a very high memory cost.

## The matrix scan

Searches the matrix will I am still filling the matrix. It
  prints out a cigar for any alignment path that is above 
  the "-min-score " setting (default 1000, which is  at
  least 200 matches) and if the path is no longer extended
  on the next row. This results in the full alignment,
  rather than best alignment being printed. So, you will
  have to search the cigars for the best sections.

The output file is named prefix--matrix-scan.aln. Each line
  in the file is tab delimited and has the score,
  ending query base, ending reference base, cigar (end of
  alignment to start), start of query, start of reference
  to a file named.
  
One problem is that it will also print out alignments that
  extend the main alignment by one base, which means this
  MAKES VERY LARGE FILES and takes a lot of time. Also, the
  print function is not very well done.

## The multi-base print

This records a best score (above "-min-score") for each
  reference base and query base. However, this will only
  record a score if the alignment is not extended into the
  next row. The alternative alignments are printed out with
  the same file format as the matrix scan.

Though this will not make as large files as the matrix
  scan, it still can output files that are the size of the
  query length * reference length. This can result in VERY
  LARGE FILES for large alignments.

This was an intresting idea, but for the few tests that I 
  ran, I found that this mainly feature printed alignments
  that were just one or two bases off the main alignment.

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

# Thanks

- To my dad for being willing to listen and give advice.
- To the people (never met) who coded baba
  [https://baba.sourceforge.net/](https://baba.sourceforge.net/).
  It gave me a great visual on how the Needleman Wunsch
  algorithm worked. Their were many other sources, but this
  was the one that was the most useful to me.
