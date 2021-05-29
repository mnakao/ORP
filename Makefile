ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas -mavx2
  LDFLAGS=-lm
  OMP_FLAGS=-qopenmp
else ifeq ($(ENV), fugaku)
  CC=fccpx
  CFLAGS=-O3 -Nclang
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  CFLAGS=-O3 -march=native -std=gnu99 -Wno-unknown-pragmas
  LDFLAGS=-lm
  OMP_FLAGS=-fopenmp
endif

serial: liborp.a
threads: liborp_threads.a
all: serial threads

liborp.a: aspl.o utils.o
	$(AR) r $@ $^
liborp_threads.a: aspl_threads.o utils_threads.o
	$(AR) r $@ $^

aspl.o: aspl.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
aspl_threads.o: aspl.c common.h parameter.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
utils.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
utils_threads.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<

clean:
	rm -f *.o *.a
