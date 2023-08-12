PREFIX=/usr/local/bin

CC=cc

CFLAGS=\
  -Wall\
  --std=c99 \
  -static\
  -O3

SOURCE=\
  alnSeqDefaults.c \
  cStrToNumberFun.c \
  seqStruct.c\
  scoresST.c\
  alnSetStruct.c\
  generalAlnFun.c\
  alnMatrixStruct.c\
  alnStruct.c\
  waterman.c\
  needleman.c\
  hirschberg.c\
  alnSeq.c
 
# Build findCoInfct
all:
	$(CC) $(CFLAGS) $(SOURCE) -o alnSeq || gcc $(CFLAGS) $(SOURCE) -o alnSeq || egcc $(CFLAGS) $(SOURCE) -o alnSeq || cc $(CFLAGS) $(SOURCE) -o alnSeq

noseqcnvt:
	$(CC) -DNOSEQCNVT $(CFLAGS) $(SOURCE) -o alnSeq

debug:
	$(CC) -Wall -DNOSEQCNVT -static --std=c99 -O0 -ggdb $(SOURCE) -o alnSeqDebug
	# Used to use -g, but -ggdb provides more info for gdb
	bash debug.sh

egcc:
	egcc $(CFLAGS) $(SOURCE) -o alnSeq
egccTwo:
	egcc $(CFLAGS) -DHIRSCHTWOBIT $(SOURCE) -o alnSeq
gcc:
	gcc $(CFLAGS) $(SOURCE) -o alnSeq
cc:
	cc $(CFLAGS) $(SOURCE) -o alnSeq

benchmark:
	$(CC) -static -Ofast $(SOURCE) -o alnSeqOFast || gcc -static -Ofast $(SOURCE) -o alnSeqOFast || egcc -static -Ofast $(SOURCE) -o alnSeqOFast || cc -static -Ofast $(SOURCE) -o alnSeqOFast
	$(CC) -static -O3    $(SOURCE) -o alnSeqO3 || gcc -static -O3 $(SOURCE) -o alnSeqO3 || egcc -static -O3 $(SOURCE) -o alnSeqO3 || cc -static -O3 $(SOURCE) -o alnSeqO3
	$(CC) -static -O2    $(SOURCE) -o alnSeqO2 || gcc -static -O2 $(SOURCE) -o alnSeqO2 || egcc -static -O2 $(SOURCE) -o alnSeqO2 || cc -static -O2 $(SOURCE) -o alnSeqO2
	$(CC) -static -O0    $(SOURCE) -o alnSeqO0 || gcc -static -O0 $(SOURCE) -o alnSeqO0 || egcc -static -O0 $(SOURCE) -o alnSeqO0 || cc -static -O0 $(SOURCE) -o alnSeqO0


clean:
	rm alnSeqDebug || printf ""; # Only thing to clean up
	rm alnSeqOFast || printf ""; # clean up benchmarking programs
	rm alnSeqO3 || printf "";
	rm alnSeqO2 || printf "";
	rm alnSeqO0 || printf "";
	rm alnSeq*.core || printf "";
	# || printf ""; is so it does not error out

install:
	mv alnSeq $(PREFIX) || printf "Unable to install alnSeq at %s\n Change this with make PREFIX=/path/to/install install\n" $(PREFIX) && exit;
	chmod a+x $(PREFIX)/alnSeq;
