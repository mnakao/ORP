#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "orp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)

static void print_help(char *argv)
{
  ERROR("%s -H hosts -S switches -R radix [-G symmetries] [-s seed] [-E]\n", argv);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix,
                     int *symmetries, int *seed, bool *assign_evenly)
{
  if(argc < 7) print_help(*argv);
  int result;
  while((result = getopt(argc,argv,"H:S:R:G:s:E"))!=-1){
    switch(result){
    case 'H':
      *hosts = atoi(optarg);
      if(*hosts < 3)
        ERROR("-H value >= 3\n");
      break;
    case 'S':
      *switches = atoi(optarg);
      if(*switches < 2)
        ERROR("-S value >= 2\n");
      break;
    case 'R':
      *radix = atoi(optarg);
      if(*radix < 3)
        ERROR("-R value >= 3\n");
      break;
    case 'G':
      *symmetries = atoi(optarg);
      if(*symmetries < 1)
        ERROR("-G value >= 1\n");
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 'E':
        *assign_evenly = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char **argv)
{
  bool assign_evenly = false;
  int hosts = -1, switches = -1, radix = -1, symmetries = 1, diameter = -1, seed = 0, lines = -1;
  long sum = -1;
  double ASPL = -1;

  set_args(argc, argv, &hosts, &switches, &radix, &symmetries, &seed, &assign_evenly);
  if(hosts%symmetries != 0){
    fprintf(stderr, "Hosts (%d) must be divisible by Symmetries (%d)\n", hosts, symmetries);
    exit(1);
  }
  else if(switches%symmetries != 0){
    fprintf(stderr, "Switches (%d) must be divisible by Symmetries (%d)\n", switches, symmetries);
    exit(1);
  }

  ORP_Srand(seed);
  int *s_degree  = malloc(sizeof(int) * switches/symmetries);
  int *h_degree  = malloc(sizeof(int) * switches/symmetries);
  int (*edge)[2] = ORP_Generate_random_s(hosts, switches, radix, assign_evenly, symmetries, &lines, h_degree, s_degree);
  printf("%d %d %d\n", hosts, switches, radix);
  ORP_Print_edge(lines, edge);
  
  return 0;
}
