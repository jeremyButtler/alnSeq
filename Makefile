PREFIX=/usr/local/bin

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
  waterman.c\
  needleman.c\
  alnSeq.c
 
# Build findCoInfct
all:
	$(CC) $(CFLAGS) $(SOURCE) -o alnSeq || gcc $(CFLAGS) $(SOURCE) -o alnSeq || egcc $(CFLAGS) $(SOURCE) -o alnSeq || cc $(CFLAGS) $(SOURCE) -o alnSeq

debug:
	$(CC) -Wall -O0 -ggdb $(SOURCE) -o alnSeqDebug
	# Used to use -g, but -ggdb provides more info for gdb
	bash debug.sh

debugOpenBsd:
	egcc -Wall -O0 -ggdb $(SOURCE) -o alnSeqDebug
	# Used to use -g, but -ggdb provides more info for gdb
	bash debug.sh

egcc:
	egcc $(CFLAGS) $(SOURCE) -o alnSeq
gcc:
	gcc $(CFLAGS) $(SOURCE) -o alnSeq
cc:
	cc $(CFLAGS) $(SOURCE) -o alnSeq

clean:
	rm alnSeqDebug; # Only thing to clean up

install:
	mv alnSeq $(PREFIX) || printf "Unable to install alnSeq at %s\n Change this with make PREFIX=/path/to/install install\n" $(PREFIX) && exit;
	chmod a+x $(PREFIX)/alnSeq;
