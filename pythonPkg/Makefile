CC=cc
CFLAGS=\
   -DDELINSSNP

python:
	CFLAGS="$(CFLAGS)" CC=$(CC) python3 setup.py install;
	rm -r build;
	rm -r dist;
	rm -r alnSeq.*egg*;

pythonlocal:
	CFLAGS="$(CFLAGS)" CC=$(CC) python3 setup.py build;
	pip3 install --user .;
	rm -r build;
	rm -r alnSeq.*egg*;
	#rm alnSeq.*.so;

# In case you cannot install with pip
libfile:
	CFLAGS="$(CFLAGS)" CC=$(CC) python3 setup.py build;
	mv build/*lib*/*.so ./;
	rm -r build;

debug:
	CFLAGS="$(CFLAGS) -O0" CC=$(CC) python3 setup-debug.py build;
	pip install --user .;
	rm -r build;
	rm -r alnSeq.*egg*;
	egdb -ex "r alnSeqPython.py" -ex "b pyAlnSeq.c:360" --args python3;   
	#gdb -ex "r alnSeqPython.py" -ex "b pyAlnSeq.c:360" --args python3;   

