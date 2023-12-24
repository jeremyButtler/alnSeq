PREFIX=/usr/local/bin

CC=cc

# Flags I want the user to be able to overwrite
CFLAGS=\
   -DDELINSSNP

# Flags I do not want the user to overwrite
CFLAGS+=\
  -Wall\
  --std=c89\
  -static\
  -O3\
  -Wno-unused-function

# -Wno-unused-function is to supress the warnings for
# some static functions I have

# These are here for the user to overwrite
#CFLAGS=-DBLANK
DEBUGFLAGS=\
   -Wall\
   --std=c89\
   -static\
   -O0\
   -ggdb\
   -Wno-unused-function\
   -DDELINSSNP\
   -DNOSEQCNVT

all:
	$(CC) $(CFLAGS) alnSeq.c -o alnSeq

python:
	CC=$(CC) make -C pythonPkg/ python;
pythonlocal:
	CC=$(CC) make -C pythonPkg/ pythonlocal;

debug:
	$(CC) $(DEBUGFLAGS) alnSeq.c -o debugAlnSeq.o
	gdb -x debugCMDs.txt debugAlnSeq.o
	# edit debugCMDs.txt to change the gdb commands

# These settings are here for quick compiling

clean:
	rm alnSeqDebug || printf ""; # Only thing to clean up
	rm alnSeqalnSeqHirschByte || printf "";
	rm alnSeqHirschTwoBit || printf "";
	rm alnSeqOFast || printf "";
	rm alnSeqO3 || printf "";
	rm alnSeqO2 || printf "";
	rm alnSeqO0 || printf "";
	rm alnSeq*.core || printf "";
	rm alnSeqByte || printf "";
	rm alnSeqTwoBit || printf "";
	rm alnSeqFast || printf "";
	rm alnSeqMid || printf "";
	# || printf ""; is so it does not error out

install:
	mv alnSeq $(PREFIX) || printf "Unable to install alnSeq at %s\n Change this with make PREFIX=/path/to/install install\n" $(PREFIX) && exit;
	chmod a+x $(PREFIX)/alnSeq;
