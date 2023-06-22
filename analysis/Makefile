all: ssw bio-alignment seq-align

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
	make -C ../benchmark;
	mv ../alnSeqOFast ./;
	mv ../alnSeqO3 ./;
	mv ../alnSeqO2 ./;
	mv ../alnSeqO0 ./;

clean:
	rm \
		-r \
		-f \
		seq-align \
		bio-alignment\
		alnSeqOFast\
		alnSeqO3\
		alnSeqO2\
		alnSeqO0\
		Complete-Striped-Smith-Waterman-Library;

# Requires emboss be installed separately
run:
	bash benchAln.sh;
	Rscript alnSeq-benchmarking.r;
