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
GRAPHS=h32s26r4.edges h80s40r6.edges h128s16r24.edges h432s120r12.edges h1024s400r10.edges h1024s1024r5.edges h1281s139r21.edges h3800s300r30.edges h4608s250r36.edges h8208s330r48.edges h10000s200r100.edges h10000s5000r10.edges h65536s3000r64.edges
serial:  simple.x simple_s.x random.x random_s.x sa.x sa_s.x hc.x sa_ref1.x
threads: threads_simple.x threads_simple_s.x threads_random.x threads_random_s.x threads_sa.x threads_sa_s.x threads_hc.x threads_sa_ref1.x
all: serial threads
simple.x: simple.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
simple_s.x: simple_s.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
random.x: random.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
random_s.x: random_s.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
sa.x: sa.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
sa_ref1.x: sa_ref1.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
sa_s.x: sa_s.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
hc.x: hc.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
threads_simple.x: simple.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_simple_s.x: simple_s.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_random.x: random.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_random_s.x: random_s.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_sa.x: sa.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_sa_ref1.x: sa_ref1.c ./liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_sa_s.x: sa_s.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_hc.x: hc.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)

test: simple.x
	for i in $(GRAPHS); do ./simple.x ./graph/$$i; done > $(TMP)
	diff $(TMP) ./graph/ASPL.txt
	rm -f $(TMP)
clean:
	rm -f *.x *.o *~
