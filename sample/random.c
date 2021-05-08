#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "orp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)

static void print_help(char *argv)
{
  ERROR("%s -H hosts -S switches -R radix [-s seed] [-E]\n", argv);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix,
                     int *seed, bool *assign_evenly)
{
  if(argc < 7) print_help(*argv);
  int result;
  while((result = getopt(argc,argv,"H:S:R:s:E"))!=-1){
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
  int hosts = -1, switches = -1, radix = -1, diameter = -1, seed = 0, lines = -1;
  long sum = -1;
  double ASPL = -1;

  set_args(argc, argv, &hosts, &switches, &radix, &seed, &assign_evenly);
  ORP_Srand(seed);
  
  int *s_degree  = malloc(sizeof(int) * switches);
  int *h_degree  = malloc(sizeof(int) * switches);
  int (*edge)[2] = ORP_Generate_random(hosts, switches, radix, assign_evenly, &lines, h_degree, s_degree);

  //  int (*adjacency)[switches] = malloc(sizeof(int) * switches * radix);
  //  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);
  //  ORP_Print_adjacency(hosts, switches, radix, s_degree, adjacency);
  printf("%d %d %d\n", hosts, switches, radix);
  ORP_Print_edge(lines, edge);
  // ORP_Print_switch(hosts, lines, edge);

  return 0;
}
