from alnSeq import alnSeqWater
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

alnList = alnSeqWater(
   ref = refSeqStr,
   query = qrySeqStr,
   noGapBool = 1,
   fullAln = 0  # not needed (prints aligned regions only)
);

print(alnList[2]); # score for alignment

print(refHeadStr); # reference header
print(alnList[0]); # aligned reference sequence

print(qryHeadStr); # query header
print(alnList[1]); # Aligned query sequence

