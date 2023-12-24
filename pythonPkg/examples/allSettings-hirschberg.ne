from alnSeq import alnSeqHirsch
import sys

# Get files and open files
refFileStr = "../../testGenomes/Small-ref.fasta";
qryFileStr = "../../testGenomes/Small-query.fasta";
refFILE = open(refFileStr, "r");
qryFILE = open(qryFileStr, "r");

refHeadStr = "";
qryHeadStr = "";
refSeqStr = "";
qrySeqStr = "";
tmpStr = "";

refHeadStr = refFILE.readline();  # Move past the header
refHeadStr = refHeadStr.replace("\n", "");
refSeqStr = ''.join(refFILE.readlines()).replace("\n", "");
   # readlines reads in all remaining lines in the file
   # ''.join converts the list to a string
      # the '' is needed to fool python into thinking it
      # is a string
   # .replace("\n", "") is to remove new lines

#re.sub("\n", "", refSeqStr); # for regular expressions

qryHeadStr = qryFILE.readline();  # move past the header
qryHeadStr = refHeadStr.replace("\n", "");
qrySeqStr = ''.join(qryFILE.readlines()).replace("\n", "");
   # readlines reads in all remaining lines in the file
   # ''.join converts the list to a string
      # the '' is needed to fool python into thinking it
      # is a string
   # .replace("\n", "") is to remove new lines

alnList = alnSeqHirsch(
   ref = refSeqStr,   # reference sequence to align
   query = qrySeqStr, # query sequence to align
   gapOpen = -10,     # gap opening penalty
   gapExtend = -1,    # gap extension penalty
   noGapBool = 1,     # Do no use gap extension penalty
   scoreMatrix = "scoring-matrix.txt" # scoring matrix
);

print(alnList[2]); # Score for the alignment

print(refHeadStr); # reference sequence header
print(alnList[0]); # Aligned reference sequence

print(qryHeadStr); # query sequence header
print(alnList[1]); # aligned query sequence

