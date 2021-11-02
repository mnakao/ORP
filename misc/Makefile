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
serial: detect_temp.x detect_switches.x histgram_s.x histgram_h.x opt_complete_graph.x add_vertex.x hmetis.x metis.x simgrid_dijkstra.x simgrid_dijkstra_dfs.x simgrid_dijkstra_random.x simgrid_wf.x simgrid_wf_dfs.x simgrid_wf_random.x
all: serial
simgrid_dijkstra.x: simgrid_dijkstra.c
	$(CC) $(CFLAGS) -o $@ $<
simgrid_dijkstra_dfs.x: simgrid_dijkstra_dfs.c
	$(CC) $(CFLAGS) -o $@ $<
simgrid_dijkstra_random.x: simgrid_dijkstra_random.c
	$(CC) $(CFLAGS) -o $@ $<
simgrid_wf.x: simgrid_wf.c
	$(CC) $(CFLAGS) -o $@ $<
simgrid_wf_dfs.x: simgrid_wf_dfs.c
	$(CC) $(CFLAGS) -o $@ $<
simgrid_wf_random.x: simgrid_wf_random.c
	$(CC) $(CFLAGS) -o $@ $<
metis.x: metis.c
	$(CC) $(CFLAGS) -o $@ $<
hmetis.x: hmetis.c
	$(CC) $(CFLAGS) -o $@ $<
add_vertex.x: add_vertex.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
opt_complete_graph.x: opt_complete_graph.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
histgram_s.x: histgram_s.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
histgram_h.x: histgram_h.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
detect_temp.x: detect_temp.c ../liborp.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
threads_detect_temp.x: detect_temp.c ../liborp_threads.a ../include/orp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) $(OMP_FLAGS)
detect_switches.x: detect_switches.c
	$(CC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -f *.x *.o *~
