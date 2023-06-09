CC=gcc

CFLAGS=\
  -Wall\
  -static\
  -O3

SOURCE=\
  cStrToNumberFun.c \
  twoBitArrays.c \
  seqStruct.c\
  scoresST.c\
  alnSetStruct.c\
  generalAlnFun.c\
  alnMatrixStruct.c\
  alnStruct.c\
  sequenceFun.c \
  waterman.c\
  needleman.c\
  alnSeq.c
 
# Build findCoInfct
all:
	$(CC) $(CFLAGS) $(SOURCE) -o alnSeq || gcc $(CFLAGS) $(SOURCE) -o alnSeq || egcc $(CFLAGS) $(SOURCE) -o alnSeq || cc $(CFLAGS) $(SOURCE) -o alnSeq

debug:
	$(CC) -Wall -O0 -ggdb $(SOURCE) -o alnSeqDebug
	# Used to use -g, but -ggdb provides more info for gdb
debugOpenBsd:
	egcc -Wall -O0 -ggdb $(SOURCE) -o alnSeqDebug
	# Used to use -g, but -ggdb provides more info for gdb
