ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -mavx2 -I../include
  OMP_FLAGS=-qopenmp
else ifeq ($(ENV), fugaku)
  CC=fccpx
  CFLAGS=-Kfast -I../include -Nclang
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  CFLAGS=-O3 -std=gnu99 -I../include
  OMP_FLAGS=-fopenmp
endif

LIBS=../liborp.a -lm
THREADS_LIBS=../liborp_threads.a -lm
TMP := $(shell mktemp)
serial:  sa_ref1.x 
threads: threads_sa_ref1.x
all: serial threads
sa_ref1.x: sa_ref1.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
threads_sa_ref1.x: sa_ref1.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)

clean:
	rm -f *.x *.o *~
