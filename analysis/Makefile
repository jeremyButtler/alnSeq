all: ssw bio-alignment alnSeq
20230722: ssw bio-alignment seq-align 20230722alnSeq

ssw:
	git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library;
	cd Complete-Striped-Smith-Waterman-Library/src;
	make -C Complete-Striped-Smith-Waterman-Library/src || gmake -C Complete-Striped-Smith-Waterman-Library/src;

bio-alignment:
	git clone https://github.com/masyagin1998/bio-alignment;
	make -C bio-alignment || gmake -C bio-alignment;

seq-align:
	git clone --recursive https://github.com/noporpoise/seq-align;
	make -C seq-align || gmake -C seq-align;

alnSeq:
	make -C ../ benchmark;
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
	rm	alnSeqO3 || printf ""\
	rm	alnSeqO2 || printf ""\
	rm	alnSeqO0 || printf ""\
	rm	alnSeqHirschByte || printf ""\
	rm	alnSeqHirschTwoBit || printf ""\

# Requires emboss be installed separately
run:
	bash benchAln.sh;
	Rscript alnSeq-benchmarking.r;
