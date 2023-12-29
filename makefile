# Makefile of STRsensor (2020/11/28)

CC = gcc
CFLAGS = -std=c99 -fopenmp
LIBS = -lz -lm -Lhtslib
HTSLIB = htslib/libhts.a
INCLUDE = -Ihtslib

DEBUG = 0

ifeq ($(DEBUG), 1)
    CFLAGS += -g -O0 # enable debugging
else
    CFLAGS += -O3
endif


OBJECT = bamio.o getseq.o index.o kmer.o cigar.o model.o parse.o utils.o main.o
PROG = STRsensor


$(PROG): $(OBJECT) $(HTSLIB)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)


# generate object file (*.o) for each source file (*.c)
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<


.PHONY : clean
clean:
	rm -f $(OBJECT)

.PHONY : fullclean
fullclean:
	rm -f $(OBJECT) $(PROG) && cd htslib && $(MAKE) clean && $(MAKE) distclean

