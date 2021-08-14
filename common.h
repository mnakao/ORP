#ifndef COMMON_INCLUDED
#define COMMON_INCLUDED

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include "parameter.h"
#if defined(__ARM_NEON)
#include <arm_neon.h>
#elif !defined(__FUJITSU)
#include <immintrin.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
  int u[2], v[2], u_d[2], v_d[2], switches, symmetries, op;
} ORP_Restore;

#define DISCONNECTED_GRAPH 0
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define UINT64_BITS       64
#define ASPL_MATRIX        1
#define ASPL_BFS           2
#define NOT_VISITED        0
#define VISITED            1
#define NOT_DEFINED       -1
#define NOT_USED          -1 
#define OP_SWAP            0
#define	OP_SWING           1

#endif
