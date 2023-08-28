all: ssw bio-alignment alnSeq
20230722: ssw bio-alignment seq-align 20230722alnSeq
20230811: ssw bio-alignment 20230811alnSeq

CC=gcc
CXX=g++

ssw:
	git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library || printf "ssw already downloaded\n";
	cd Complete-Striped-Smith-Waterman-Library/src || printf "error happened when downloading ssw\n";
	make CC=$(CC) CXX=$(CXX) -C Complete-Striped-Smith-Waterman-Library/src || gmake CC=$(CC) CXX=$(CXX) -C Complete-Striped-Smith-Waterman-Library/src;

bio-alignment:
	git clone https://github.com/masyagin1998/bio-alignment;
	make CC=$(CC) -C bio-alignment || gmake CC=$(CC) -C bio-alignment;

seq-align:
	git clone --recursive https://github.com/noporpoise/seq-align;
	make CC=$(CC) -C seq-align || gmake CC=$(CC) -C seq-align;

alnSeq:
	make -C ../ benchmark;
	mv ../alnSeqTwoBit ./;
	mv ../alnSeqByte ./;
	mv ../alnSeqMid ./;
	mv ../alnSeqFast ./;

20230811alnSeq:
	make -C ../ 20230812benchmark;
	mv ../alnSeqHirschByte ./;
	mv ../alnSeqHirschTwoBit ./;


20230811alnSeq:
	make -C ../ 20230812benchmark;
	mv ../alnSeqHirschByte ./;
	mv ../alnSeqHirschTwoBit ./;

20230722alnSeq:
	make -C ../ 20230722benchmark;
	mv ../alnSeqOFast ./;
	mv ../alnSeqO3 ./;
	mv ../alnSeqO2 ./;
	mv ../alnSeqO0 ./;

clean:
	rm -r -f bio-alignment || printf "";
	rm -r -f Complete-Striped-Smith-Waterman-Library || printf "";
	rm -r -f seq-align || printf "";
	rm alnSeqOFast || printf "";
	rm	alnSeqO3 || printf "";
	rm	alnSeqO2 || printf "";
	rm	alnSeqO0 || printf "";
	rm	alnSeqHirschByte || printf "";
	rm	alnSeqHirschTwoBit || printf "";
	rm alnSeqByte || printf "";
	rm alnSeqTwoBit || printf "";
	rm alnSeqFast || printf "";
	rm alnSeqMid || printf "";

# Requires emboss be installed separately
run:
	bash 20230827-benchmarking/20230827-benchAln.sh;
	Rscript 20230827-benchmarking/20230827-alnSeq-benchmarking.r;

run20230811:
	bash 20230811-benchmarking/20230811-benchAln.sh;
	Rscript 20230811-benchmarking/20230811-alnSeq-benchmarking.r;

run20230722:
	bash 20230722-benchmarking/20230722-benchAln.sh;
	Rscript 20230722-benchmarking/20230722-alnSeq-benchmarking.r;
