#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "orp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define NOT_DEFINED -1
#define DEFAULT_SEED 0
#define DEFAULT_NCALCS 100

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

bool accept(const int switches, const int current_diameter, const int diameter, const double current_ASPL,
            const double ASPL, double *max_diff_energy)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }

  //  diameter == current_diameter
  double diff = (current_ASPL-ASPL)*switches*switches;
  *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);
    
  return (ASPL <= current_ASPL);
}

static void print_help(char *argv)
{
  fprintf(stderr, "%s [-H hosts] [-S switches] [-R radix] [-s seed] [-A] [-E]\n", argv);
  fprintf(stderr, "  -H : Number of hosts (Required when -f option is not specified)\n");
  fprintf(stderr, "  -S : Number of switches (Set automatically when -f and -S are not specified).\n");
  fprintf(stderr, "  -R : Radix (Required when -f is not specified)\n");
  fprintf(stderr, "  -s : Random seed (Default: %d)\n", DEFAULT_SEED);
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix, int *seed)
{
  int result;
  while((result = getopt(argc,argv,"H:S:R:s:"))!=-1){
    switch(result){
    case 'H':
      *hosts = atoi(optarg);
      if(*hosts <= 3)
        ERROR("-H value > 3\n");
      break;
    case 'S':
      *switches = atoi(optarg);
      if(*switches <= 2)
        ERROR("-S value > 2\n");
      break;
    case 'R':
      *radix = atoi(optarg);
      if(*radix <= 3)
        ERROR("-R value > 3\n");
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, seed = DEFAULT_SEED;
  int lines, diameter, current_diameter, low_diameter;
  int (*edge)[2], *h_degree, *s_degree;
  long sum, best_sum, ncalcs = DEFAULT_NCALCS;
  double ASPL, current_ASPL, low_ASPL, max_diff_energy = 0;
  ORP_Restore r;

  set_args(argc, argv, &hosts, &switches, &radix, &seed);
  
  ORP_Srand(seed);
  if(hosts == NOT_DEFINED || radix == NOT_DEFINED)
    print_help(argv[0]);
  else if(switches == NOT_DEFINED)
    switches = ORP_Optimize_switches(hosts, radix);
  
  h_degree = malloc(sizeof(int) * switches);
  s_degree = malloc(sizeof(int) * switches);
  edge     = ORP_Generate_random(hosts, switches, radix, false, &lines, h_degree, s_degree);
  
  printf("Hosts = %d, Switches = %d, Radix = %d\n", hosts, switches, radix);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);

  int (*adjacency)[radix]      = malloc(sizeof(int) * switches * radix);

  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);
  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  current_diameter = diameter;
  current_ASPL     = ASPL;

  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    for(long i=0;i<ncalcs;i++){
      if(uniform_rand() > 0.5)
        ORP_Swap_adjacency(switches, radix, s_degree, &r, adjacency);
      else
        ORP_Swing_adjacency(switches, radix, h_degree, s_degree, &r, adjacency);
      
      ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
      
      if(accept(switches, current_diameter, diameter, current_ASPL, ASPL, &max_diff_energy)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
        ORP_Restore_adjacency(r, radix, h_degree, s_degree, adjacency);
      }
    }
  }
  ORP_Finalize_aspl();

  // Set max temperature to accept it   50% in maximum diff energy.
  printf("Proposed max temperature is %f\n", (-1.0 * max_diff_energy) / log(0.5));
  // Set min temperature to accept it 0.01% in minimum diff energy.
  printf("Proposed min temperature is %f\n", (-2.0) / log(0.0001));
  
  return 0;
}
