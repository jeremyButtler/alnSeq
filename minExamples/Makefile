CC=cc
CFLAGS=\
   -static\
   -DDELINSSNP\
   -Wno-unsued-function\
   -O3

all:
	$(CC) $(CFLAGS) -o needle-no-gap-extend.o needle-no-gap-extend.c
clean:
	rm needle-no-gap-extend.o || printf "";