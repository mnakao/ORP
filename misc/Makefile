ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -mavx2 -I../include
  OMP_FLAGS=-qopenmp
else ifeq ($(ENV), fugaku)
  CC=fccpx
  CFLAGS=-O3 -I../include -Nclang
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  CFLAGS=-O3 -std=gnu99 -I../include
  OMP_FLAGS=-fopenmp
endif

LIBS=../liborp.a -lm
THREADS_LIBS=../liborp_threads.a -lm
TMP := $(shell mktemp)
serial: detect_temp.x detect_switches.x histgram.x
all: serial
histgram.x: histgram.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
detect_temp.x: detect_temp.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
threads_detect_temp.x: detect_temp.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) $(OMP_FLAGS)
detect_switches.x: detect_switches.c
	$(CC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -f *.x *.o *~
